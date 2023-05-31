

// "Wind Walker Herd" by dr2 - 2023
// License: Creative Commons Attribution-NonCommercial-ShareAlike 4.0

#define AA  0   // (= 0/1) optional antialiasing

#if 0
#define VAR_ZERO min (nFrame, 0)
#else
#define VAR_ZERO 0
#endif

uniform float iTime;
float iFrame = iTime;
uniform vec2 iResolution;
uniform vec2 iMouse;

float PrBoxDf (vec3 p, vec3 b);
float PrBox2Df (vec2 p, vec2 b);
float PrRoundBoxDf (vec3 p, vec3 b, float r);
float PrCylDf (vec3 p, float r, float h);
float PrCylAnDf (vec3 p, float r, float w, float h);
float PrCaps2Df (vec2 p, float r, float h);
float Minv2 (vec2 p);
float Maxv2 (vec2 p);
mat3 StdVuMat (float el, float az);
vec2 Rot2D (vec2 q, float a);
vec2 Rot2Cs (vec2 q, vec2 cs);
float ShowIntPZ (vec2 q, vec2 cBox, float mxChar, float val);
vec2 Hashv2v2 (vec2 p);
float Fbm1 (float p);
float Fbm2 (vec2 p);
vec3 VaryNf (vec3 p, vec3 n, float f);

vec3 sunDir, qHit, lBase;
vec2 cId;
float tCur, dstFar, bGrid, grLim, spd, wlkScl;
int nFrame, idObj;
const int idBas = 1, idLnkT = 2, idLnk = 3, idAx = 4, idWhl = 5, idVane = 6, idStruc = 7,
   idGrnd = 8;
const float s_a = 38.0, s_b = 41.5, s_c = 39.3, s_d = 40.1, s_e = 55.8, s_f = 39.4,
   s_g = 36.7, s_h = 65.7, s_i = 49.0, s_j = 50.0, s_k = 61.9, s_l = 7.8, s_m = 15.0;
const float pi = 3.1415927;

#define CosSin(x) (sin ((x) + vec2 (0.5 * pi, 0.)))
#define DMINQ(id) if (d < dMin) { dMin = d;  idObj = id;  qHit = q; }

struct Leg {
  vec2 v[8], cs[10], cswAng;
  float wAng;
};
struct Walker {
  Leg leg[2];
  vec2 csVane;
  float szFac;
};
Walker wlk;

#define ACOSR(x, y, z) acos (((x) * (x) + (y) * (y) - (z) * (z)) / (2. * (x) * (y)))
#define ATANV(v) atan ((v).y, (v).x)

void ObjState ()
{ //  (Leg from the Strandbeest: see https://en.wikipedia.org/wiki/Jansen's_linkage)
  float a[10], aa, g, s, t;
  wlk.szFac = wlkScl / (1. + 0.3 * Maxv2 (abs (cId)));
  t = tCur * wlkScl / wlk.szFac;
  wlk.leg[0].wAng = - spd * t;
  wlk.leg[1].wAng = wlk.leg[0].wAng + pi;
  for (int k = 0; k < 2; k ++) {
    wlk.leg[k].v[0] = vec2 (0., 0.);
    wlk.leg[k].v[1] = wlk.leg[k].v[0] + vec2 (s_a, s_l);
    wlk.leg[k].v[2] = wlk.leg[k].v[1] + Rot2D (vec2 (s_m, 0.), wlk.leg[k].wAng);
    aa = ATANV (wlk.leg[k].v[2] - wlk.leg[k].v[0]);
    s = length (wlk.leg[k].v[2] - wlk.leg[k].v[0]);
    a[0] = aa + ACOSR (s, s_b, s_j);
    wlk.leg[k].v[3] = wlk.leg[k].v[0] + Rot2D (vec2 (s_b, 0.), a[0]);
    a[1] = aa - ACOSR (s, s_c, s_k);
    wlk.leg[k].v[4] = wlk.leg[k].v[0] + Rot2D (vec2 (s_c, 0.), a[1]);
    a[2] = ACOSR (s_b, s_d, s_e) + a[0];
    wlk.leg[k].v[5] = wlk.leg[k].v[0] + Rot2D (vec2 (s_d, 0.), a[2]);
    s = length (wlk.leg[k].v[5] - wlk.leg[k].v[4]);
    g = ACOSR (s, s_c, s_d) + ACOSR (s, s_g, s_f) + pi + a[1];
    wlk.leg[k].v[6] = wlk.leg[k].v[4] + Rot2D (vec2 (s_g, 0.), g);
    wlk.leg[k].v[7] = wlk.leg[k].v[4] + Rot2D (vec2 (s_i, 0.), g + ACOSR (s_g, s_i, s_h));
    a[3] = ACOSR (s_d, s_e, s_b) + a[2] - pi;
    a[4] = ATANV (wlk.leg[k].v[4] - wlk.leg[k].v[6]);
    a[5] = ATANV (wlk.leg[k].v[5] - wlk.leg[k].v[6]);
    a[6] = ATANV (wlk.leg[k].v[7] - wlk.leg[k].v[6]);
    a[7] = ATANV (wlk.leg[k].v[7] - wlk.leg[k].v[4]);
    a[8] = ATANV (wlk.leg[k].v[3] - wlk.leg[k].v[2]);
    a[9] = ATANV (wlk.leg[k].v[4] - wlk.leg[k].v[2]);
    for (int m = 0; m < 10; m ++) wlk.leg[k].cs[m] = CosSin (- a[m]);
    wlk.leg[k].cswAng = CosSin (- wlk.leg[k].wAng);
  }
  wlk.csVane = CosSin (4. * t);
}

void LinkDf (vec3 p, vec2 v, vec2 cs, float l, int id, inout float dMin)
{
  vec3 q;
  float d;
  q = p;
  q.xy = Rot2Cs (q.xy - v, cs);
  d = max (PrCaps2Df (q.yx - vec2 (0., 0.5 * l), 2., 0.5 * l), abs (q.z) - 0.5);
  DMINQ (id);
}

float ObjDf (vec3 p)
{ // (Based on "Wind Walker")
  vec3 q, pp;
  float dMin, d, sx;
  p.xz -= bGrid * (cId + 0.5);
  dMin = dstFar / wlk.szFac;
  p /= wlk.szFac;
  p.y -= lBase.y;
  p.xz = Rot2Cs (p.xz, CosSin (0.25 * pi));
  p.xz = p.zx * vec2 (1., -1.);
  pp = p;
  for (int kx = 0; kx < 2; kx ++) {
    sx = sign (float (kx) - 0.5);
    p.x = pp.x  + lBase.x * sx;
    for (int k = 0; k < 2; k ++) {
      p.z = lBase.z + pp.z * (sign (float (k) - 0.5)) * sx;
      q = p;
      q.z -= 4.;
      q.xy = Rot2Cs (q.xy - wlk.leg[k].v[1], wlk.leg[k].cswAng);
      d = PrCylAnDf (q, s_m, 2., 1.);
      q.xy = (abs (q.x) > abs (q.y)) ? q.xy : q.yx;
      d = min (d, max (PrBox2Df (q.xy, vec2 (s_m, 1.8)), abs (q.z) - 0.8));
      DMINQ (idWhl);
      for (int j = 2; j <= 7; j ++) {
        q = p;
        q -= vec3 (wlk.leg[k].v[j], 0.2);
        d = PrCylDf (q, ((j < 7) ? 1.5 : 2.5), ((j == 2) ? 5. : 3.));
        DMINQ (idAx);
      }
      LinkDf (p, wlk.leg[k].v[0], wlk.leg[k].cs[0], s_b, idLnkT, dMin);
      LinkDf (p, wlk.leg[k].v[0], wlk.leg[k].cs[2], s_d, idLnkT, dMin);
      LinkDf (p, wlk.leg[k].v[5], wlk.leg[k].cs[3], s_e, idLnkT, dMin);
      LinkDf (p, wlk.leg[k].v[6], wlk.leg[k].cs[4], s_g, idLnkT, dMin);
      LinkDf (p, wlk.leg[k].v[6], wlk.leg[k].cs[6], s_h, idLnkT, dMin);
      LinkDf (p, wlk.leg[k].v[4], wlk.leg[k].cs[7], s_i, idLnkT, dMin);
      p.z -= 1.4;
      LinkDf (p, wlk.leg[k].v[0], wlk.leg[k].cs[1], s_c, idLnk, dMin);
      LinkDf (p, wlk.leg[k].v[6], wlk.leg[k].cs[5], s_f, idLnk, dMin);
      LinkDf (p, wlk.leg[k].v[2], wlk.leg[k].cs[8], s_j, idLnk, dMin);
      p.z += 2.8;
      LinkDf (p, wlk.leg[k].v[2], wlk.leg[k].cs[9], s_k, idLnk, dMin);
    }
  }
  p = pp;
  q = p;
  q.x -= 20.;
  d = PrRoundBoxDf (q, vec3 (lBase.x + 35., 2.5, lBase.z - 7.), 0.5);
  DMINQ (idBas);
  q = p;
  q.x = abs (q.x) - lBase.x;
  d = PrCylDf (q, 1.5, lBase.z + 2.);
  DMINQ (idAx);
  q = p;
  q.xy -= vec2 (s_a, s_l);
  q.x = abs (q.x) - lBase.x;
  d = PrCylDf (q, 1.5, lBase.z - 2.);
  DMINQ (idAx);
  q = p;
  q.xy -= vec2 (s_a, s_l - 1.5);
  d = PrRoundBoxDf (q, vec3 (lBase.x, 4.5, 4.), 0.5);
  DMINQ (idStruc);
  q.x = abs (abs (q.x) - 0.5 * lBase.x) - 0.5 * lBase.x;
  d = PrRoundBoxDf (q, vec3 (4., 4.5, lBase.z - 7.), 0.5);
  DMINQ (idStruc);
  q = p;
  q.z = abs (abs (q.z) - 21.);
  q -= vec3 (s_a, 27., 21.);
  d = PrBoxDf (q, vec3 (4., 27., 1.5));
  DMINQ (idStruc);
  q = p;
  q.xy -= vec2 (s_a, 50.);
  d = PrCylDf (q, 2.5, lBase.z - 5.);
  DMINQ (idAx);
  q.xy = Rot2Cs (q.xy, wlk.csVane);
  d = max (abs (length (q.xy - vec2 (18., 10.)) - 20.) - 0.2, q.y);
  q.xy = Rot2Cs (q.xy, CosSin (2. * pi / 3.));
  d = min (d, max (abs (length (q.xy - vec2 (18., 10.)) - 20.) - 0.2, q.y));
  q.xy = Rot2Cs (q.xy, CosSin (2. * pi / 3.));
  d = min (d, max (abs (length (q.xy - vec2 (18., 10.)) - 20.) - 0.2, q.y));
  q.z = abs (q.z) - 21.;
  d = max (d, abs (q.z) - 18.);
  DMINQ (idVane);
  return wlk.szFac * dMin;
}

float ObjRay (vec3 ro, vec3 rd)
{
  vec3 p, rdi;
  vec2 cIdP, s;
  float dHit, d, eps;
  eps = 0.0005;
  if (rd.x == 0.) rd.x = 0.001;
  if (rd.z == 0.) rd.z = 0.001;
  rdi.xz = 1. / rd.xz;
  cIdP = vec2 (-999.);
  dHit = eps;
  for (int j = VAR_ZERO; j < 160; j ++) {
    p = ro + dHit * rd;
    cId = floor (p.xz / bGrid);
    if (cId != cIdP) {
      ObjState ();
      cIdP = cId;
    }
    d = (Maxv2 (abs (cId)) <= grLim) ? ObjDf (p) : dstFar;
    s = (bGrid * (cId + step (0., rd.xz)) - p.xz) * rdi.xz;
    d = min (d, abs (Minv2 (s)) + eps);
    dHit += d;
    if (d < eps || dHit > dstFar || p.y < 0.) break;
  }
  if (d >= eps) dHit = dstFar;
  return dHit;
}

vec3 ObjNf (vec3 p)
{
  vec4 v;
  vec2 e;
  e = vec2 (0.001, -0.001);
  for (int j = VAR_ZERO; j < 4; j ++) {
    v[j] = ObjDf (p + ((j < 2) ? ((j == 0) ? e.xxx : e.xyy) : ((j == 2) ? e.yxy : e.yyx)));
  }
  v.x = - v.x;
  return normalize (2. * v.yzw - dot (v, vec4 (1.)));
}

float ObjSShadow (vec3 ro, vec3 rd)
{
  vec3 p;
  vec2 cIdP;
  float sh, d, h;
  sh = 1.;
  d = 0.02;
  cIdP = vec2 (-999.);
  for (int j = VAR_ZERO; j < 24; j ++) {
    p = ro + d * rd;
    cId = floor (p.xz / bGrid);
    if (cId != cIdP) {
      ObjState ();
      cIdP = cId;
    }
    if (Maxv2 (abs (cId)) <= grLim) {
      h = ObjDf (p);
      sh = min (sh, smoothstep (0., 0.05 * d, h));
    } else h = 0.3 * bGrid;
    d += clamp (h, 0.02, 0.5);
    if (sh < 0.05 || d > dstFar) break;
  }
  return 0.6 + 0.4 * sh;
}

vec3 SkyBgCol (vec3 ro, vec3 rd)
{
  vec3 col, clCol, skCol;
  vec2 q;
  float f, fd, ff, sd;
  if (rd.y > -0.02 && rd.y < 0.03 * Fbm1 (16. * atan (rd.z, - rd.x))) {
    col = vec3 (0.3, 0.41, 0.55);
  } else if (rd.y < 0.) {
    col = vec3 (0.3, 0.41, 0.55);
  } else {
    q = 0.02 * (ro.xz + 0.5 * tCur + ((100. - ro.y) / rd.y) * rd.xz);
    ff = Fbm2 (q);
    f = smoothstep (0.2, 0.8, ff);
    fd = smoothstep (0.2, 0.8, Fbm2 (q + 0.01 * sunDir.xz)) - f;
    clCol = (0.7 + 0.5 * ff) * (vec3 (0.7) - 0.7 * vec3 (0.3, 0.3, 0.2) * sign (fd) *
       smoothstep (0., 0.05, abs (fd)));
    sd = max (dot (rd, sunDir), 0.);
    skCol = vec3 (0.4, 0.5, 0.8) + step (0.1, sd) * vec3 (1., 1., 0.9) *
       min (0.3 * pow (sd, 64.) + 0.5 * pow (sd, 2048.), 1.);
    col = mix (skCol, clCol, 0.1 + 0.9 * f * smoothstep (0.01, 0.1, rd.y));
  }
  return 0.8 * col;
}

vec4 ObjCol ()
{
  vec4 col4;
  if (idObj == idBas) {
    col4 = vec4 (0.8, 0.6, 0.2, 0.05) * (0.9 +
       0.1 * smoothstep (0.1, 0.13, fract (8. * abs (qHit.z) / 50. + 0.5)));
    if (qHit.y > 0.) col4 = mix (col4, vec4 (0.2, 1., 0.2, -1.), ShowIntPZ (Rot2D (qHit.xz +
       vec2 (20. + lBase.x, 20.), 0.5 * pi), 0.7 * lBase.x * vec2 (1., 0.5),
       2., 11. + grLim - cId.x + (2. * grLim + 1.) * (grLim - cId.y)));
  } else if (idObj == idStruc) {
    col4 = vec4 (0.5, 0.5, 0.8, 0.05);
  } else if (idObj == idLnkT) {
    col4 = vec4 (0.85, 0.85, 0.9, 0.1) * (0.8 + 0.2 * smoothstep (0.18, 0.22, abs (qHit.y)));
  } else if (idObj == idLnk) {
    col4 = vec4 (0.95, 0.95, 1., 0.1) * (0.8 + 0.2 * smoothstep (0.18, 0.22,
       abs (abs (qHit.y) - 0.8)));
  } else if (idObj == idAx) {
    col4 = vec4 (0.8, 0.7, 0.2, 0.1);
  } else if (idObj == idWhl) {
    col4 = vec4 (1., 0.6, 0.2, 0.05);
  } else if (idObj == idVane) {
    col4 = vec4 (1., 1., 0.9, 0.05) * (0.8 + 0.2 * smoothstep (0.2, 0.24,
       abs (abs (qHit.z) - 10.)));
  }
  return col4;
}

vec3 ShowScene (vec3 ro, vec3 rd)
{
  vec4 col4;
  vec3 col, vn, q;
  float dstObj, dstGrnd, sh, t, nDotL;
  bool isBg;
  wlkScl = 0.07;
  lBase = vec3 (60., 86., 50.);
  spd = 2.;
  dstGrnd = dstFar;
  isBg = false;
  dstObj = ObjRay (ro, rd);
  if (dstObj < dstFar) {
    ro += dstObj * rd;
    vn = ObjNf (ro);
    col4 = ObjCol ();
  } else if (rd.y < 0.) {
    dstGrnd = - ro.y / rd.y;
    ro += dstGrnd * rd;
    q = ro;
    q.xz += spd * tCur;
    col4 = mix (vec4 (0.4, 0.5, 0.3, 0.), vec4 (0., 0.5, 0.1, 0.), smoothstep (0.2, 0.8, Fbm2 (q.xz)));
    col4 = mix (vec4 (0.2, 0.5, 0.2, 0.), col4,  1. - smoothstep (0.6, 0.9, dstGrnd / dstFar));
    vn = VaryNf (2. * q, vec3 (0., 1., 0.), 2. * (1. - smoothstep (0.2, 0.4, dstGrnd / dstFar)));
  } else {
    col = SkyBgCol (ro, rd);
    isBg = true;
  }
  if (! isBg) {
    if (col4.a >= 0.) {
      nDotL = max (dot (vn, sunDir), 0.);
      if (dstObj < dstFar && (idObj == idLnk || idObj == idLnkT || idObj == idWhl)) nDotL *= nDotL;
      sh = (min (dstObj, dstGrnd) < dstFar) ? ObjSShadow (ro + 0.01 * vn, sunDir) : 1.;
      col = col4.rgb * (0.2 + 0.2 * max (dot (vn, sunDir * vec3 (-1., 1., -1.)), 0.) +
         0.8 * sh * nDotL) + step (0.95, sh) * col4.a * pow (max (0.,
         dot (sunDir, reflect (rd, vn))), 32.);
    } else col = col4.rgb * (0.55 - 0.45 * dot (rd, vn));
    if (dstObj >= dstFar) col = mix (col, 0.8 * vec3 (0.3, 0.41, 0.55), pow (1. + rd.y, 16.));
  }
  return clamp (col, 0., 1.);
}

void main() {
  mat3 vuMat;
  vec4 mPtr;
  vec3 ro, rd, col;
  vec2 canvas, uv;
  float el, az, zmFac, sr;
  nFrame = iFrame;
  canvas = iResolution.xy;
  uv = 1. * gl_FragCoord.xy / canvas - 1.;
  uv.x *= canvas.x / canvas.y;
  tCur = iTime;
  mPtr.xy = iMouse;
  mPtr.xy = mPtr.xy / canvas - 0.5;
  bGrid = 16.;
  grLim = 3.;
  az = 0.;
  el = -0.1 * pi;

// Uncomment following lines to enable click to drag
//   if (mPtr.z > 0.) {
    az -= 2. * pi * mPtr.x;
    el -= 0.5 * pi * mPtr.y;
//   } else {
//     az = mod (az + 0.01 * pi * tCur + pi, 2. * pi) - pi;
//   }
  el = clamp (el, -0.4 * pi, -0.03 * pi);
  vuMat = StdVuMat (el, az);
  ro = vuMat * vec3 (0., 0., -150.);
  ro.xz += 0.5 * bGrid;
  zmFac = 6. + 3. * abs (az);
  dstFar = 300.;
  sunDir = normalize (vec3 (0., 1., -1.));
  sunDir.xz = Rot2D (sunDir.xz, -0.01 * pi * tCur);
#if ! AA
  const float naa = 1.;
#else
  const float naa = 3.;
#endif
  col = vec3 (0.);
  sr = 2. * mod (dot (mod (floor (0.5 * (uv + 1.) * canvas), 2.), vec2 (1.)), 2.) - 1.;
  for (float a = float (VAR_ZERO); a < naa; a ++) {
    rd = vuMat * normalize (vec3 (uv + step (1.5, naa) * Rot2D (vec2 (0.5 / canvas.y, 0.),
       sr * (0.667 * a + 0.5) * pi), zmFac));
    col += (1. / naa) * ShowScene (ro, rd);
  }
  gl_FragColor = vec4 (col, 1.);
}

float PrBoxDf (vec3 p, vec3 b)
{
  vec3 d;
  d = abs (p) - b;
  return min (max (d.x, max (d.y, d.z)), 0.) + length (max (d, 0.));
}

float PrBox2Df (vec2 p, vec2 b)
{
  vec2 d;
  d = abs (p) - b;
  return min (max (d.x, d.y), 0.) + length (max (d, 0.));
}

float PrRoundBoxDf (vec3 p, vec3 b, float r)
{
  return length (max (abs (p) - b, 0.)) - r;
}

float PrCylDf (vec3 p, float r, float h)
{
  return max (length (p.xy) - r, abs (p.z) - h);
}

float PrCylAnDf (vec3 p, float r, float w, float h)
{
  return max (abs (length (p.xy) - r) - w, abs (p.z) - h);
}

float PrCaps2Df (vec2 p, float r, float h)
{
  return length (vec2 (p.x, sign (p.y) * (max (0., abs (p.y) - h)))) - r;
}

float Minv2 (vec2 p)
{
  return min (p.x, p.y);
}

float Maxv2 (vec2 p)
{
  return max (p.x, p.y);
}

mat3 StdVuMat (float el, float az)
{
  vec2 ori, ca, sa;
  ori = vec2 (el, az);
  ca = cos (ori);
  sa = sin (ori);
  return mat3 (ca.y, 0., - sa.y, 0., 1., 0., sa.y, 0., ca.y) *
         mat3 (1., 0., 0., 0., ca.x, - sa.x, 0., sa.x, ca.x);
}

vec2 Rot2D (vec2 q, float a)
{
  vec2 cs;
  cs = sin (a + vec2 (0.5 * pi, 0.));
  return vec2 (dot (q, vec2 (cs.x, - cs.y)), dot (q.yx, cs));
}

vec2 Rot2Cs (vec2 q, vec2 cs)
{
  return vec2 (dot (q, vec2 (cs.x, - cs.y)), dot (q.yx, cs));
}

float DigSeg (vec2 q)
{
  q = 1. - smoothstep (vec2 (0.), vec2 (0.04, 0.07), abs (q) - vec2 (0.13, 0.5));
  return q.x * q.y;
}

#define DSG(q) k = kk;  kk = k / 2;  if (kk * 2 != k) d += DigSeg (q)

float ShowDig (vec2 q, int iv)
{
  vec2 vp, vm, vo;
  float d;
  int k, kk;
  vp = vec2 (0.5, 0.5);
  vm = vec2 (-0.5, 0.5);
  vo = vp - vm;
  if (iv == -1) k = 8;
  else if (iv < 2) k = (iv == 0) ? 119 : 36;
  else if (iv < 4) k = (iv == 2) ? 93 : 109;
  else if (iv < 6) k = (iv == 4) ? 46 : 107;
  else if (iv < 8) k = (iv == 6) ? 122 : 37;
  else             k = (iv == 8) ? 127 : 47;
  q = (q - 0.5) * vec2 (1.8, 2.3);
  d = 0.;
  kk = k;
  DSG (q.yx - vo);  DSG (q.xy - vp);  DSG (q.xy - vm);  DSG (q.yx);
  DSG (q.xy + vm);  DSG (q.xy + vp);  DSG (q.yx + vo);
  return d;
}

float ShowIntPZ (vec2 q, vec2 cBox, float mxChar, float val)
{
  float nDig, idChar, s, v;
  q = vec2 (- q.x, q.y) / cBox;
  s = 0.;
  if (Minv2 (q) >= 0. && Maxv2 (q) < 1.) {
    q.x *= mxChar;
    nDig = mxChar;
    idChar = mxChar - 1. - floor (q.x);
    q.x = fract (q.x);
    v = max (val, 0.) / pow (10., mxChar - idChar - 1.);
    if (idChar >= mxChar - nDig) s = ShowDig (q, int (mod (floor (v), 10.)));
  }
  return s;
}

const float cHashM = 43758.54;

vec2 Hashv2f (float p)
{
  return fract (sin (p + vec2 (0., 1.)) * cHashM);
}

vec2 Hashv2v2 (vec2 p)
{
  vec2 cHashVA2 = vec2 (37., 39.);
  return fract (sin (vec2 (dot (p, cHashVA2), dot (p + vec2 (1., 0.), cHashVA2))) * cHashM);
}

float Noiseff (float p)
{
  vec2 t;
  float ip, fp;
  ip = floor (p);
  fp = fract (p);
  fp = fp * fp * (3. - 2. * fp);
  t = Hashv2f (ip);
  return mix (t.x, t.y, fp);
}

float Noisefv2 (vec2 p)
{
  vec2 t, ip, fp;
  ip = floor (p);
  fp = fract (p);
  fp = fp * fp * (3. - 2. * fp);
  t = mix (Hashv2v2 (ip), Hashv2v2 (ip + vec2 (0., 1.)), fp.y);
  return mix (t.x, t.y, fp.x);
}

float Fbm1 (float p)
{
  float f, a;
  f = 0.;
  a = 1.;
  for (int j = 0; j < 5; j ++) {
    f += a * Noiseff (p);
    a *= 0.5;
    p *= 2.;
  }
  return f * (1. / 1.9375);
}

float Fbm2 (vec2 p)
{
  float f, a;
  f = 0.;
  a = 1.;
  for (int j = 0; j < 5; j ++) {
    f += a * Noisefv2 (p);
    a *= 0.5;
    p *= 2.;
  }
  return f * (1. / 1.9375);
}

float Fbmn (vec3 p, vec3 n)
{
  vec3 s;
  float a;
  s = vec3 (0.);
  a = 1.;
  for (int j = 0; j < 5; j ++) {
    s += a * vec3 (Noisefv2 (p.yz), Noisefv2 (p.zx), Noisefv2 (p.xy));
    a *= 0.5;
    p *= 2.;
  }
  return dot (s, abs (n));
}

vec3 VaryNf (vec3 p, vec3 n, float f)
{
  vec3 g;
  vec2 e = vec2 (0.1, 0.);
  g = vec3 (Fbmn (p + e.xyy, n), Fbmn (p + e.yxy, n), Fbmn (p + e.yyx, n)) - Fbmn (p, n);
  return normalize (n + f * (g - n * dot (n, g)));
}
