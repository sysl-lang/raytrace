import struct, sys, zlib
import numpy as np

# Decode an 8-bit RGB or RGBA PNG with zlib and numpy alone, undoing the five row filters, and print
# each channel's mean and standard deviation.
data = open(sys.argv[1], "rb").read()
pos, idat, w, h, ct = 8, b"", 0, 0, 0
while pos < len(data):
    n, t = struct.unpack(">I4s", data[pos:pos + 8])
    body = data[pos + 8:pos + 8 + n]
    if t == b"IHDR":
        w, h, depth, ct = struct.unpack(">IIBB", body[:10])
    elif t == b"IDAT":
        idat += body
    pos += 12 + n
ch = {2: 3, 6: 4}[ct]
raw = zlib.decompress(idat)
stride = w * ch
out = np.zeros((h, stride), dtype=np.int64)
prev = np.zeros(stride, dtype=np.int64)
for y in range(h):
    f = raw[y * (stride + 1)]
    line = np.frombuffer(raw[y * (stride + 1) + 1:(y + 1) * (stride + 1)], dtype=np.uint8).astype(np.int64)
    cur = np.zeros(stride, dtype=np.int64)
    for x in range(stride):
        a = cur[x - ch] if x >= ch else 0
        b = prev[x]
        c = prev[x - ch] if x >= ch else 0
        if f == 0: p = 0
        elif f == 1: p = a
        elif f == 2: p = b
        elif f == 3: p = (a + b) // 2
        else:
            pa, pb, pc = abs(b - c), abs(a - c), abs(a + b - 2 * c)
            p = a if pa <= pb and pa <= pc else (b if pb <= pc else c)
        cur[x] = (line[x] + p) & 255
    out[y] = cur
    prev = cur
img = out.reshape(h, w, ch)[:, :, :3].astype(np.float64)
print(w, h)
for k in range(3):
    print("mean", repr(img[:, :, k].mean()), "sd", repr(img[:, :, k].std()))
