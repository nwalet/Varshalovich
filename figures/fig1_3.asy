// Fig. 1.3 -- succession of rotations (scheme A), three Euler stages.
// One shared orthographic viewpoint (x comes forward, matches Figs 1.1/1.2);
// a single continuous frame carried across the three panels.
//   asy -f pdf fig1_3.asy   ->  fig1_3a.pdf, fig1_3b.pdf, fig1_3c.pdf
import three;
settings.outformat="pdf"; settings.prc=false; settings.render=4;
currentprojection=orthographic((4,1.2,2.6));   // x forward, ~30 deg elevation

real R=2, L=2.6, Lz=3.0;
pen kpen=black+linewidth(0.9);
pen bpen=blue+linewidth(1.0);
pen gpen=deepgreen+linewidth(1.0);
pen dpen=gray(0.45)+linewidth(0.7)+dashed;     // carried-over (grey) axis
pen dbpen=blue+linewidth(0.7)+dashed;          // carried-over (blue) axis
pen greyfill=gray(0.85)+opacity(0.42);
pen bluefill=lightblue+opacity(0.);
pen redfill=red+opacity(0.2);

// continuous frame -------------------------------------------------------
transform3 Rz=rotate(-35,Z);            triple x1=Rz*X,  y1=Rz*Y;   // stage a (z1=Z)
transform3 Ryb=rotate(10,y1);          triple z2=Ryb*Z, x2=Ryb*x1; // stage b (y2=y1)
transform3 Rz2=rotate(30,z2);          triple xp=Rz2*x2, yp=Rz2*y1;// stage c

surface diskS(triple n){ return surface(circle(O,R,n)); }
void fixbox(picture pic){                // identical bounding box -> identical scale/framing
  for(int i=-1;i<=1;i+=2) for(int j=-1;j<=1;j+=2) for(int k=-1;k<=1;k+=2)
    dot(pic,(3.2*i,3.2*j,3.2*k),invisible);
}
void vaxis(picture pic, triple v, real len, pen p, string s, align al){
  draw(pic, O--len*v, p, Arrow3(size=5)); label(pic, s, len*v, al, black);
}
void refaxis(picture pic, triple v, real len, pen p, string s, align al, pen lc){
  draw(pic, O--len*v, p); label(pic, s, len*v, al, lc);
}
void anglearc(picture pic, triple va, triple vb, real rad, triple ax, string s,
              align al=N, real lf=1.35){
  draw(pic, arc(O, rad*unit(va), rad*unit(vb), ax), black+linewidth(0.7), Arrow3(size=4));
  label(pic, s, lf*rad*unit(va+vb), al, black);
}

// ---------- panel (a): rotate about z by alpha ----------
picture pa; fixbox(pa);
draw(pa, diskS(Z), greyfill,grey);
vaxis(pa, X, L, kpen, "$x$", S);
vaxis(pa, Y, L, kpen, "$y$", E);
vaxis(pa, Z, Lz, kpen, "$z=z_1$", N);
vaxis(pa, x1, L, bpen, "$x_1$", S);
vaxis(pa, y1, L, bpen, "$y_1$", N);
anglearc(pa, X, x1, 1.15, -Z, "$\alpha$", S, 1.5);

// ---------- panel (b): rotate about y1 by beta ----------
picture pb; fixbox(pb);
draw(pb, diskS(Z),  white+opacity(0.0), grey);                   // original xy-plane, as in (a)
draw(pb, diskS(z2), white+opacity(1.0));                    // tilted plane
draw(pb, diskS(y1), redfill, red);                    // tilted plane
refaxis(pb, Z,  Lz, dpen, "$z_1$", N, gray(0.4));
refaxis(pb, x1, L,  dpen, "$x_1$", S, gray(0.4));
vaxis(pb, y1, L,  kpen, "$y_1=y_2$", NE);
vaxis(pb, z2, Lz, bpen, "$z_2$", NE);
vaxis(pb, x2, L,  bpen, "$x_2$", S);
anglearc(pb, Z, z2, 1.35, y1, "$\beta$", W, 1.15);

// ---------- panel (c): rotate about z2=z' by gamma ----------
picture pc; fixbox(pc);
draw(pc, diskS(Z),  greyfill);
draw(pc, diskS(z2), bluefill);
refaxis(pc, x2, L, dbpen, "$x_2$", S, blue);
refaxis(pc, y1, L, dbpen, "$y_2$", E, blue);
vaxis(pc, z2, Lz, kpen, "$z_2=z'$", NW);
vaxis(pc, xp, L, gpen, "$x'$", S);
draw(pc, O--L*yp, gpen, Arrow3(size=5)); label(pc, "$y'$", 1.16*L*yp, E, black);
anglearc(pc, x2, xp, 1.35, z2, "$\gamma$", E, 1.2);

shipout("fig1_3a", pa.fit(5cm));
shipout("fig1_3b", pb.fit(5cm));
shipout("fig1_3c", pc.fit(5cm));
