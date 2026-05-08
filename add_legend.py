from PIL import Image, ImageDraw, ImageFont
import colorsys

Image.MAX_IMAGE_PIXELS = None

def hex_color(h, s, l):
    r, g, b = colorsys.hls_to_rgb(h, l, s)
    return (int(r*255), int(g*255), int(b*255))

def shade(tidx, h, s):
    l = 0.88 - (min(tidx, 14) / 14.0) * 0.45
    return hex_color(h, s, l)

H_TAG=0.58; H_SEC=0.38; H_U=0.07; H_HYST=0.97; H_PRED=0.72

block_types = [
    (H_TAG,  0.55, "Tag"),
    (H_SEC,  0.50, "Sec Tag"),
    (H_U,    0.75, "U"),
    (H_HYST, 0.65, "Hys"),
    (H_PRED, 0.60, "Pred"),
]

SPECIAL = [
    ((0, 168, 150), "FB Pred"),
    ((0, 120, 110), "FB Hys"),
    ((244, 196, 48), "Meta"),
]

img = Image.open("floorplan_final.png").convert("RGBA")
W, H_img = img.size

scale = max(W, H_img) / 1800.0

MARGIN    = int(18  * scale)
PADDING   = int(14  * scale)
ROW_H     = int(32  * scale)
SWATCH_W  = int(28  * scale)
SWATCH_H  = int(20  * scale)
FONT_SIZE = int(20  * scale)
GRAD_H    = int(22  * scale)
GRAD_W    = int(260 * scale)
GAP       = int(10  * scale)
SHIFT     = int(0  * scale)

try:
    font    = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", FONT_SIZE)
    font_sm = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", max(FONT_SIZE - 3, 8))
    font_hd = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", FONT_SIZE)
except:
    font = font_sm = font_hd = ImageFont.load_default()

tmp = ImageDraw.Draw(Image.new("RGBA", (1, 1)))
max_label_w = max(int(tmp.textlength(label, font=font)) for _, _, label in block_types)
max_label_w = max(max_label_w, max(int(tmp.textlength(label, font=font)) for _, label in SPECIAL))

legend_w = max(PADDING*2 + SWATCH_W + GAP + max_label_w, GRAD_W + PADDING*2) + int(30 * scale)
n_rows   = len(block_types) + len(SPECIAL)
legend_h = (PADDING*2
            + FONT_SIZE + GAP
            + n_rows * ROW_H + GAP
            + FONT_SIZE + GAP
            + GRAD_H + GAP
            + FONT_SIZE)

lx = W     - MARGIN - legend_w - SHIFT
ly = H_img - MARGIN - legend_h

overlay = Image.new("RGBA", img.size, (0, 0, 0, 0))
draw = ImageDraw.Draw(overlay)
draw.rounded_rectangle(
    [lx, ly, lx + legend_w, ly + legend_h],
    radius=int(12 * scale),
    fill=(15, 15, 15, 210)
)
img = Image.alpha_composite(img, overlay)
draw = ImageDraw.Draw(img)

y = ly + PADDING

draw.text((lx + PADDING, y), "Block type", font=font_hd, fill=(255, 255, 255, 255))
y += FONT_SIZE + GAP

for h, s, label in block_types:
    color = shade(7, h, s)
    sx = lx + PADDING
    sy = y + (ROW_H - SWATCH_H) // 2
    draw.rectangle([sx, sy, sx + SWATCH_W, sy + SWATCH_H], fill=color + (255,))
    draw.text((sx + SWATCH_W + GAP, y + (ROW_H - FONT_SIZE) // 2), label, font=font, fill=(235, 235, 235, 255))
    y += ROW_H

for color, label in SPECIAL:
    sx = lx + PADDING
    sy = y + (ROW_H - SWATCH_H) // 2
    draw.rectangle([sx, sy, sx + SWATCH_W, sy + SWATCH_H], fill=color + (255,))
    draw.text((sx + SWATCH_W + GAP, y + (ROW_H - FONT_SIZE) // 2), label, font=font, fill=(235, 235, 235, 255))
    y += ROW_H

y += GAP
draw.text((lx + PADDING, y), "Table index  t0 → t14", font=font_sm, fill=(200, 200, 200, 255))
y += FONT_SIZE + GAP

step_w = GRAD_W // 15
for i in range(15):
    color = shade(i, H_TAG, 0.55)
    x0 = lx + PADDING + i * step_w
    draw.rectangle([x0, y, x0 + step_w, y + GRAD_H], fill=color + (255,))

y += GRAD_H + GAP
draw.text((lx + PADDING, y), "t0", font=font_sm, fill=(180, 180, 180, 255))
t14_w = int(tmp.textlength("t14", font=font_sm))
draw.text((lx + PADDING + GRAD_W - t14_w, y), "t14", font=font_sm, fill=(180, 180, 180, 255))

img = img.convert("RGB")
img.save("floorplan_legend.png")
print("Saved floorplan_legend.png")
