import math, struct, sys, zlib

out = sys.argv[1]

# An icosahedron subdivided once and pushed onto the unit sphere: 42 vertices, 80 faces.
t = (1 + 5 ** 0.5) / 2
v = [(-1, t, 0), (1, t, 0), (-1, -t, 0), (1, -t, 0), (0, -1, t), (0, 1, t), (0, -1, -t), (0, 1, -t),
     (t, 0, -1), (t, 0, 1), (-t, 0, -1), (-t, 0, 1)]
v = [tuple(c / math.sqrt(sum(x * x for x in p)) for c in p) for p in v]
f = [(0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11), (1, 5, 9), (5, 11, 4), (11, 10, 2),
     (10, 7, 6), (7, 1, 8), (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9), (4, 9, 5), (2, 4, 11),
     (6, 2, 10), (8, 6, 7), (9, 8, 1)]
mid = {}
def m(a, b):
    k = (min(a, b), max(a, b))
    if k not in mid:
        p = [(v[a][i] + v[b][i]) / 2 for i in range(3)]
        n = math.sqrt(sum(x * x for x in p))
        v.append(tuple(c / n for c in p))
        mid[k] = len(v) - 1
    return mid[k]
nf = []
for a, b, c in f:
    ab, bc, ca = m(a, b), m(b, c), m(c, a)
    nf += [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
with open(out + "/icosphere.obj", "w") as o:
    o.write("# An icosahedron subdivided once onto the unit sphere: 42 vertices, 80 triangles.\n")
    for p in v:
        o.write("v %.6f %.6f %.6f\n" % p)
    for a, b, c in nf:
        o.write("f %d %d %d\n" % (a + 1, b + 1, c + 1))

# A 128x64 planet: latitude bands of ocean and land, a grid every 30 degrees, polar caps.
W, H = 128, 64
rows = []
for y in range(H):
    row = bytearray([0])
    lat = 90 - (y + 0.5) * 180 / H
    for x in range(W):
        lon = (x + 0.5) * 360 / W - 180
        land = math.sin(math.radians(lon) * 3) * math.cos(math.radians(lat) * 2) > 0.2
        r, g, b = (60, 140, 50) if land else (20, 60, 170)
        if abs(lat) > 70:
            r, g, b = 235, 240, 245
        if abs(lon % 30) < 1.5 or abs(lat % 30) < 1.5:
            r, g, b = 250, 220, 80
        row += bytes([r, g, b])
    rows.append(bytes(row))
def chunk(t, d):
    return struct.pack(">I", len(d)) + t + d + struct.pack(">I", zlib.crc32(t + d) & 0xffffffff)
png = b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", W, H, 8, 2, 0, 0, 0)) + \
    chunk(b"IDAT", zlib.compress(b"".join(rows), 9)) + chunk(b"IEND", b"")
open(out + "/planet.png", "wb").write(png)
