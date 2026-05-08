import re, colorsys, sys

def hex_color(h, s, l):
    r, g, b = colorsys.hls_to_rgb(h, l, s)
    return "#{:02x}{:02x}{:02x}".format(int(r*255), int(g*255), int(b*255))

def shade(tidx, h, s):
    l = 0.88 - (min(tidx, 14) / 14.0) * 0.45
    return hex_color(h, s, l)

H_TAG=0.58; H_SEC=0.38; H_U=0.07; H_HYST=0.97; H_PRED=0.72
SPECIAL_FB="#00a896"; SPECIAL_FB_BIM="#00786e"; SPECIAL_META="#f4c430"

input_file  = sys.argv[1] if len(sys.argv) > 1 else "floorplan.gv"
output_file = sys.argv[2] if len(sys.argv) > 2 else "floorplan_final.gv"

with open(input_file) as f:
    lines = f.readlines()

rect_info = {}
anchors = []
for line in lines:
    m = re.match(r'\s*rect(\d+)\s*\[label="\d+(?:\\n([^"]*))?"', line)
    if m:
        rnum = int(m.group(1))
        name = m.group(2) or ""
        rect_info[rnum] = name
        tm = re.match(r't(\d+)_(\w+)', name)
        if tm:
            anchors.append((rnum, int(tm.group(1)), tm.group(2)))

def get_color(rnum):
    name = rect_info.get(rnum, "")
    if 'fb_bim_hyst' in name: return SPECIAL_FB_BIM
    if name.startswith('fb'):  return SPECIAL_FB
    if name.startswith('meta'): return SPECIAL_META
    tm = re.match(r't(\d+)_(\w+)', name)
    if tm:
        tidx, btype = int(tm.group(1)), tm.group(2)
    else:
        best = min(anchors, key=lambda a: abs(a[0] - rnum))
        tidx, btype = best[1], best[2]
    colors = {"tag":(H_TAG,0.55),"sec":(H_SEC,0.50),
              "u":(H_U,0.75),"hyst":(H_HYST,0.65),"pred":(H_PRED,0.60)}
    h, s = colors.get(btype, (0.6, 0.4))
    return shade(tidx, h, s)

out = []
for line in lines:
    # Keep only name, drop the number
    line = re.sub(r'(label=")(\d+)\\n([^"]*)"', lambda m: m.group(1)+m.group(3)+'"', line)
    line = re.sub(r'(label=")\d+"', r'\1"', line)
    # Replace fillcolor
    m = re.match(r'(\s*rect(\d+)\s*\[.*?fillcolor=")[^"]*(".*)', line)
    if m:
        line = m.group(1) + get_color(int(m.group(2))) + m.group(3) + "\n"
    out.append(line)

with open(output_file, "w") as f:
    f.writelines(out)
print(f"Written {output_file}")
