#include <iostream>
#include <fstream>
#include <iomanip>
#include "trace_reader.hpp"

static constexpr uint64_t MAX_INSTRUCTIONS = 100000;

const char* inst_class_name(INST_CLASS c) {
    switch (c) {
        case INST_CLASS::ALU:                return "ALU";
        case INST_CLASS::LOAD:               return "LOAD";
        case INST_CLASS::STORE:              return "STORE";
        case INST_CLASS::BR_COND:            return "BR_COND";
        case INST_CLASS::BR_UNCOND_DIRECT:   return "BR_UNCOND_DIRECT";
        case INST_CLASS::BR_UNCOND_INDIRECT: return "BR_UNCOND_INDIRECT";
        case INST_CLASS::FP:                 return "FP";
        case INST_CLASS::ALU_SLOW:           return "ALU_SLOW";
        case INST_CLASS::BR_CALL_DIRECT:     return "BR_CALL_DIRECT";
        case INST_CLASS::BR_CALL_INDIRECT:   return "BR_CALL_INDIRECT";
        case INST_CLASS::BR_RETURN:          return "BR_RETURN";
        default:                             return "UNKNOWN";
    }
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " trace.gz output.txt [limit]\n";
        std::cerr << "  limit: max instructions to decode (default: "
                  << MAX_INSTRUCTIONS << ")\n";
        return 1;
    }

    uint64_t limit = MAX_INSTRUCTIONS;
    if (argc > 3) limit = std::stoull(argv[3]);

    trace_reader reader(argv[1], "decode");
    std::ofstream out(argv[2]);

    if (!out.is_open()) {
        std::cerr << "Failed to open output file: " << argv[2] << "\n";
        return 1;
    }

    out << "# Trace decode: " << argv[1] << "\n";
    out << "# Limit: " << limit << " instructions\n";
    out << "#\n";
    out << "# Format: INDEX PC=0x... NEXT=0x... TYPE=... [TAKEN=Y/N]\n";
    out << "#\n";

    uint64_t count = 0;
    try {
        while (count < limit) {
            instruction inst = reader.next_instruction();

            out << std::dec << std::right << std::setw(10) << count
                << " PC=0x"   << std::hex << std::setw(16) << std::setfill('0') << inst.pc
                << " NEXT=0x" << std::setw(16) << inst.next_pc
                << std::setfill(' ') << std::dec
                << " TYPE="   << std::left << std::setw(20)
                << inst_class_name(inst.inst_class);

            if (inst.branch)
                out << " TAKEN=" << (inst.taken_branch ? "Y" : "N");

            out << "\n";
            count++;
        }
    }
    catch (const out_of_instructions&) { /* normal EOF */ }

    std::cout << "Decoded " << count << " instructions to " << argv[2] << "\n";
    if (count == limit)
        std::cout << "(limit reached — trace may have more instructions)\n";

    return 0;
}