//
// Portions copyright (c) 2025 Ampere Computing
//
// This file handles reading a subset of the format used for CBP2025 at
// https://github.com/ramisheikh/cbp2025/blob/main/lib/trace_reader.h. The
// original copyright/license message from that file is reproduced below:

// CBP Trace Reader
// Author: Arthur Perais (arthur.perais@gmail.com) for CVP
//         Saransh Jain/Rami Sheikh updated for CBP

/*
   This is free and unencumbered software released into the public domain.

   Anyone is free to copy, modify, publish, use, compile, sell, or
   distribute this software, either in source code form or as a compiled
   binary, for any purpose, commercial or non-commercial, and by any
   means.

   In jurisdictions that recognize copyright laws, the author or authors
   of this software dedicate any and all copyright interest in the
   software to the public domain. We make this dedication for the benefit
   of the public at large and to the detriment of our heirs and
   successors. We intend this dedication to be an overt act of
   relinquishment in perpetuity of all present and future rights to this
   software under copyright law.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
   IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
   OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
   ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
   OTHER DEALINGS IN THE SOFTWARE.

   For more information, please refer to <http://unlicense.org>
*/

#pragma once

#include <cassert>
#include <cstdint>
#include <zlib.h>
#include <string>
#include <stdexcept>
#include <vector>
#include <algorithm>
#include <optional>

enum class INST_CLASS: uint8_t 
{
    ALU = 0,
    LOAD = 1,
    STORE = 2,
    BR_COND = 3,
    BR_UNCOND_DIRECT = 4,
    BR_UNCOND_INDIRECT = 5,
    FP = 6,
    ALU_SLOW = 7,
    UNDEF = 8,
    BR_CALL_DIRECT = 9,
    BR_CALL_INDIRECT = 10,
    BR_RETURN = 11
};

struct instruction {
    uint64_t pc = 0xdeadc0dedeadc0de;
    uint64_t next_pc = 0xdeadc0dedeadc0de;
    INST_CLASS inst_class = INST_CLASS::ALU;
    bool branch = false;
    bool taken_branch = false;
};

struct out_of_instructions: std::runtime_error {
    using std::runtime_error::runtime_error;
};


class trace_reader
{
private:
    std::string trace_name;
    gzFile gz_fp;
    // branch_only stores whether there is information in this trace about
    // non-branch instructions.
    bool branch_only;
    std::optional<struct instruction> putback_instruction;
    uint64_t num_instructions_read = 0;

public:
    trace_reader(std::string path, std::string name): trace_name{name} {
        gz_fp = gzopen(path.c_str(), "rb");
        if (!gz_fp) {
            throw std::runtime_error("Failed to open gzipped file: " + path);
        }

        // Try to read a magic number as the first 8 bytes of this trace to see
        // which version it is.
        const char magic[9] = "CBPNGAmp";
        branch_only = true;
        for (unsigned i = 0; i < sizeof(magic) - 1; i++) {
            char c;
            read(c);
            if (c != magic[i]) {
                branch_only = false;
                gzseek(gz_fp, 0, SEEK_SET);
                break;
            }
        }
    }

    std::string name() { return trace_name; }

    ~trace_reader() {
        if (gz_fp) {
            gzclose(gz_fp);
            gz_fp = nullptr;
        }
    }

private:
    template <typename T>
    void read(T& obj) {
        long unsigned num_read = gzread(gz_fp, static_cast<void*>(&obj), sizeof(T));
        if (num_read < sizeof(T)) {
            if (num_read == 0 && gzeof(gz_fp)) {
                throw out_of_instructions("Encountered EOF when attempting to read instruction");
            }
            throw std::runtime_error("Failed to read all bytes from trace");
        }
    }

    void ignore_read(size_t bytes) {
        auto ret = gzseek(gz_fp, bytes, SEEK_CUR);
        if (ret < 0)
            throw std::runtime_error("Failed to ignore bytes from trace");
    }

    struct instruction read_next_instruction() {
        // Trace Format :
        // Inst PC                  - 8 bytes
        // Inst Type                - 1 byte
        // If load/store
        //   Effective Address      - 8 bytes
        //   Access Size (total)    - 1 byte
        //   Involves Base Update   - 1 byte
        //   If Store:
        //      Involves Reg Offset - 1 byte
        // If branch
        //   Taken                  - 1 byte
        //   If Taken:
        //      Target              - 8 bytes
        // Num Input Regs           - 1 byte
        // Input Reg Names          - 1 byte each
        // Num Output Regs          - 1 byte
        // Output Reg Names         - 1 byte each
        // Output Reg Values
        //   If INT                 - 8 bytes each
        //   If SIMD                - 16 bytes each
        //
        // Int registers are encoded 0-30(GPRs), 31(Stack Pointer Register), 64(Flag Register), 65(Zero Register)
        // SIMD registers are encoded 32-63

        struct instruction instr;

        read(instr.pc);
        instr.next_pc = instr.pc + 4;

        read(instr.inst_class);

        assert(instr.inst_class != INST_CLASS::UNDEF);

        //EffAddr is the base address
        bool has_base_update = false;
        if (instr.inst_class == INST_CLASS::LOAD || instr.inst_class == INST_CLASS::STORE) {
            //   Effective Address      - 8 bytes
            //   Access Size (total)    - 1 byte
            ignore_read(8 + 1);
            //   Involves Base Update   - 1 byte
            read(has_base_update);
            //   If Store:
            //      Involves Reg Offset - 1 byte
            if (instr.inst_class == INST_CLASS::STORE)
                ignore_read(1);
        }

        const bool is_branch = 
            instr.inst_class == INST_CLASS::BR_COND || 
            instr.inst_class == INST_CLASS::BR_UNCOND_DIRECT || 
            instr.inst_class == INST_CLASS::BR_UNCOND_INDIRECT || 
            instr.inst_class == INST_CLASS::BR_CALL_DIRECT || 
            instr.inst_class == INST_CLASS::BR_CALL_INDIRECT || 
            instr.inst_class == INST_CLASS::BR_RETURN;


        if (is_branch) {
            instr.branch = true;
            read(instr.taken_branch);
            if(instr.inst_class != INST_CLASS::BR_COND)
                assert(instr.taken_branch);
            if(instr.taken_branch)
                read(instr.next_pc);
        }

        uint8_t num_in_regs, num_out_regs;
        // note: will only hold integer registers in these vectors!
        std::vector<uint8_t> in_regs, out_regs;

        read(num_in_regs);
        for (int i = 0; i < num_in_regs; i++) {
            uint8_t reg_no;
            read(reg_no);
            if (reg_no < 32)
                in_regs.push_back(reg_no);
        }

        read(num_out_regs);
        for (int i = 0; i < num_out_regs; i++) {
            uint8_t reg_no;
            read(reg_no);
            out_regs.push_back(reg_no);
        }

        // Determine if this instruction has a "base update register" so that
        // we can correctly read the trace later depending on the answer
        uint8_t base_update_reg = 0;
        if (has_base_update) {
            if (instr.inst_class == INST_CLASS::STORE) {
                if (num_out_regs != 1) {
                    has_base_update = false;
                } else {
                    assert(out_regs.size() == 1);
                    base_update_reg = out_regs[0];
                }
            } else { // loads
                if (num_out_regs <= 1) {
                    has_base_update = false;
                } else {
                    std::vector<uint8_t> overlap;
                    std::vector<uint8_t> int_out_regs;
                    std::copy_if(out_regs.begin(), out_regs.end(), std::back_inserter(int_out_regs), [](int reg) { return reg < 32; });
                    std::set_intersection(in_regs.begin(), in_regs.end(), int_out_regs.begin(), int_out_regs.end(), std::back_inserter(overlap));
                    if (overlap.size() == 1) {
                        base_update_reg = overlap[0];
                    } else {
                        has_base_update = false;
                    }
                }
            }
        }

        // ignore reads for all register values
        for (auto i = 0; i < num_out_regs; i++)
        {
            ignore_read(8);

            // Base update registers are not written to the trace, so we
            // shouldn't ignore them
            const bool matching_base_upd = has_base_update && base_update_reg == out_regs[i];
            // < 32 are integer registers; 32-63 are FP/SIMD; 64 and 65 are condition code and zero register, respectively
            const bool integer_register = (out_regs[i] < 32) || (out_regs[i] == 64) || (out_regs[i] == 65);
            if (!matching_base_upd && !integer_register)
                ignore_read(8);
        }

        num_instructions_read++;
        return instr;
    }

public:
    void put_back(struct instruction inst) {
        assert(!putback_instruction.has_value());
        putback_instruction = inst;
    }

    struct instruction next_instruction() {
        if (putback_instruction.has_value()) {
            struct instruction tmp = putback_instruction.value();
            putback_instruction.reset();
            return tmp;
        } else {
            return read_next_instruction();
        }
    }
};
