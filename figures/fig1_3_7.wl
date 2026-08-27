(* ::Package:: *)

 R=1.5;L=1.5;


x={1,0,0};y={0,1,0};z={0,0,1};


degrad=\[Pi]/180.


{x1,y1,z1}=(r1=RotationTransform[15 degrad,z])/@{x,y,z}


{x2,y2,z2}=(r2=RotationTransform[25 degrad,y1])/@{x1,y1,z1}


{x3,y3,z3}=(r2=RotationTransform[20 degrad,z2])/@{x2,y2,z2}


pla=Graphics3D[{{Directive[White,EdgeForm[Gray]],ResourceFunction["Disk3D"][{0,0,0},1,{0,0,1}]},{Red,Arrow[{{0,0,0},R x}],Arrow[{{0,0,0},R y}],Arrow[{{0,0,0},R z}]},{Blue,Arrow[{{0,0,0},R x1}],Arrow[{{0,0,0},R y1}]},
{Blue,Arrowheads[0.03],Arrow[Tube[{x,x1},0.015]],Arrow[Tube[{y,y1},0.015]]},
{Black,Arrowheads[0.03],Arrow[BezierCurve[Table[0.8 R z+0.2 {Cos[\[Theta]],Sin[\[Theta]],0},{\[Theta],-0.4,\[Pi],0.1}]],-0.1]},Text[Style["\!\(\*
StyleBox[\"z\",\nFontSlant->\"Italic\"]\)=\!\(\*SubscriptBox[
StyleBox[\"z\",\nFontSlant->\"Italic\"], \(1\)]\)",20,FontFamily->"Times New Roman"],R z,{0,-1}],Text[Style["\!\(\*
StyleBox[\"x\",\nFontSlant->\"Italic\"]\)",20,FontFamily->"Times New Roman"],R x,{1,0}],Text[Style["\!\(\*SubscriptBox[\(x\), \(1\)]\)",20,FontFamily->"Times New Roman"],R x1,{0,1}],Text[Style["\!\(\*
StyleBox[\"y\",\nFontSlant->\"Italic\"]\)",20,FontFamily->"Times New Roman"],R y,{-1,0}],Text[Style["\!\(\*SubscriptBox[\(y\), \(1\)]\)",20,FontFamily->"Times New Roman"],R y1,{-1,0}],Text[Style["(a)",20,FontFamily->"Times New Roman"],{0.5L,-0.8L,L}]},Lighting->"Neutral",ViewPoint->{10,10,10},Boxed->False,PlotRange->L{{-1,1},{-1,1},{-1,1}}]


plb=Graphics3D[{{Directive[Opacity[0],EdgeForm[{Dashed,Black}]],ResourceFunction["Disk3D"][{0,0,0},1,{0,0,1}]},{Directive[White,Opacity[0.7],EdgeForm[{Blue}]],ResourceFunction["Disk3D"][{0,0,0},1,z2]},{Directive[Red,EdgeForm[{Red}]],ResourceFunction["Disk3D"][{0,0,0},1,y2]},{Blue,Arrow[{{0,0,0},R x1}],Arrow[{{0,0,0},R y1}],Arrow[{{0,0,0},R z1}]},{Darker[Green],Arrow[{{0,0,0},R x2}],Arrow[{{0,0,0},R z2}]},{Green,Arrowheads[0.03],Arrow[Tube[{x1,x2},0.015]],Arrow[Tube[{z1,z2},0.015]]},{Black,Arrowheads[0.03],Arrow[BezierCurve[Reverse[Table[0.8 R y1+r1[0.2 {Cos[\[Theta]],0,Sin[\[Theta]]}],{\[Theta],-0.4,\[Pi],0.1}]]],-0.1]},Text[Style["\!\(\*SubscriptBox[\(y\), \(1\)]\)=\!\(\*SubscriptBox[\(y\), \(2\)]\)",20,FontFamily->"Times New Roman"],R y1,{-1,0}],Text[Style["\!\(\*SubscriptBox[\(x\), \(1\)]\)",20,FontFamily->"Times New Roman",Background->White],R x1,{1,0}],Text[Style["\!\(\*SubscriptBox[\(x\), \(2\)]\)",20,FontFamily->"Times New Roman"],R x2,{0,1}],Text[Style["\!\(\*SubscriptBox[\(z\), \(1\)]\)",20,FontFamily->"Times New Roman"],R z1,{0,-1}],Text[Style["\!\(\*SubscriptBox[\(z\), \(2\)]\)",20,FontFamily->"Times New Roman"],R z2,{0,-1}],Text[Style["(b)",20,FontFamily->"Times New Roman"],{0.5L,-0.8L,L}]},Lighting->"Neutral",ViewPoint->{10,10,6},Boxed->False,PlotRange->L{{-1,1},{-1,1},{-1,1}}]


plc=Graphics3D[{{Directive[Opacity[0],EdgeForm[{Dashed,Black}]],ResourceFunction["Disk3D"][{0,0,0},1,{0,0,1}]},{Directive[White,Opacity[0.7],EdgeForm[{Blue}]],ResourceFunction["Disk3D"][{0,0,0},1,z2]},{Directive[Opacity[0],EdgeForm[{Red}]],ResourceFunction["Disk3D"][{0,0,0},1,y2]},{Directive[Darker[Green],Opacity[0.5],EdgeForm[{Green}]],ResourceFunction["Disk3D"][{0,0,0},1,y3]},{Darker[Green],Arrow[{{0,0,0},R x2}],Arrow[{{0,0,0},R y2}],Arrow[{{0,0,0},R z2}]},{Purple,Arrow[{{0,0,0},R x3}],Arrow[{{0,0,0},R y3}]},{Purple,Arrowheads[0.03],Arrow[Tube[{x2,x3},0.015]],Arrow[Tube[{y2,y3},0.015]]},{Black,Arrowheads[0.03],Arrow[BezierCurve[Table[0.8 R z2+r2[r1[0.2 {Cos[\[Theta]],Sin[\[Theta]],0}]],{\[Theta],-0.4,\[Pi],0.1}]],-0.1]},Text[Style[TraditionalForm["\!\(\*
StyleBox[\"z\",\nFontSlant->\"Italic\"]\)'=\!\(\*SubscriptBox[\(z\), \(2\)]\)"],20,FontFamily->"Times New Roman"],R z2,{0,-1}],Text[Style["\!\(\*SubscriptBox[\(x\), \(2\)]\)",20,FontFamily->"Times New Roman"],R x2,{1,0}],Text[Style["\!\(\*
StyleBox[\"x\",\nFontSlant->\"Italic\"]\)'",20,FontFamily->"Times New Roman"],R x3,{0,1}],Text[Style[TraditionalForm["\!\(\*SubscriptBox[\(y\), \(2\)]\)"],20,FontFamily->"Times New Roman"],R y2,{-1,0}],Text[Style[TraditionalForm["\!\(\*
StyleBox[\"y\",\nFontSlant->\"Italic\"]\)'"],20,FontFamily->"Times New Roman"],R y3,{-1,0}],Text[Style["(c)",20,FontFamily->"Times New Roman"],{0.5L,-0.8L,L}]},Lighting->"Neutral",ViewPoint->{10,10,6},Boxed->False,PlotRange->L{{-1,1},{-1,1},{-1,1}}]


Export["fig1_3.pdf",GraphicsRow[{pla,plb,plc},0,ImageSize->1024]]


{x1,y1,z1}=(r1=RotationTransform[20 degrad,z])/@{x,y,z}


{x2,y2,z2}=(r2=RotationTransform[25 degrad,y])/@{x1,y1,z1}


{x3,y3,z3}=(r2=RotationTransform[15 degrad,z])/@{x2,y2,z2}


plb=Graphics3D[{{Directive[Opacity[0],EdgeForm[{Dashed,Black}]],ResourceFunction["Disk3D"][{0,0,0},1,{0,0,1}]},{Directive[White,Opacity[0.7],EdgeForm[{Blue}]],ResourceFunction["Disk3D"][{0,0,0},1,z2]},{Directive[Red,EdgeForm[{Red}]],ResourceFunction["Disk3D"][{0,0,0},1,y2]},{Blue,Arrow[{{0,0,0},R x1}],Arrow[{{0,0,0},R y1}],Arrow[{{0,0,0},R z1}]},{Darker[Green],Arrow[{{0,0,0},R x2}],Arrow[{{0,0,0},R z2}]},{Green,Arrowheads[0.03],Arrow[Tube[{x1,x2},0.015]],Arrow[Tube[{z1,z2},0.015]]},{Black,Arrowheads[0.03],Arrow[BezierCurve[Reverse[Table[0.8 R y1+r1[0.2 {Cos[\[Theta]],0,Sin[\[Theta]]}],{\[Theta],-0.4,\[Pi],0.1}]]],-0.1]},Text[Style["\!\(\*SubscriptBox[\(y\), \(1\)]\)=\!\(\*SubscriptBox[\(y\), \(2\)]\)",20,FontFamily->"Times New Roman"],R y1,{-1,0}],Text[Style["\!\(\*SubscriptBox[\(x\), \(1\)]\)",20,FontFamily->"Times New Roman",Background->White],R x1,{1,0}],Text[Style["\!\(\*SubscriptBox[\(x\), \(2\)]\)",20,FontFamily->"Times New Roman"],R x2,{0,1}],Text[Style["\!\(\*SubscriptBox[\(z\), \(1\)]\)",20,FontFamily->"Times New Roman"],R z1,{0,-1}],Text[Style["\!\(\*SubscriptBox[\(z\), \(2\)]\)",20,FontFamily->"Times New Roman"],R z2,{0,-1}],Text[Style["(b)",20,FontFamily->"Times New Roman"],{0.5L,-0.8L,L}]},Lighting->"Neutral",ViewPoint->{10,10,6},Boxed->False,PlotRange->L{{-1,1},{-1,1},{-1,1}}]


pla=Graphics3D[{{Directive[White,EdgeForm[Gray]],ResourceFunction["Disk3D"][{0,0,0},1,{0,0,1}]},{Red,Arrow[{{0,0,0},R x}],Arrow[{{0,0,0},R y}],Arrow[{{0,0,0},R z}]},{Blue,Arrow[{{0,0,0},R x1}],Arrow[{{0,0,0},R y1}]},
{Blue,Arrowheads[0.03],Arrow[Tube[{x,x1},0.015]],Arrow[Tube[{y,y1},0.015]]},
{Black,Arrowheads[0.03],Arrow[BezierCurve[Table[0.8 R z+0.2 {Cos[\[Theta]],Sin[\[Theta]],0},{\[Theta],-0.4,\[Pi],0.1}]],-0.1]},Text[Style["\!\(\*
StyleBox[\"z\",\nFontSlant->\"Italic\"]\)=\!\(\*SubscriptBox[
StyleBox[\"z\",\nFontSlant->\"Italic\"], \(1\)]\)",20,FontFamily->"Times New Roman"],R z,{0,-1}],Text[Style["\!\(\*
StyleBox[\"x\",\nFontSlant->\"Italic\"]\)",20,FontFamily->"Times New Roman"],R x,{1,0}],Text[Style["\!\(\*SubscriptBox[\(x\), \(1\)]\)",20,FontFamily->"Times New Roman"],R x1,{0,1}],Text[Style["\!\(\*
StyleBox[\"y\",\nFontSlant->\"Italic\"]\)",20,FontFamily->"Times New Roman"],R y,{-1,0}],Text[Style["\!\(\*SubscriptBox[\(y\), \(1\)]\)",20,FontFamily->"Times New Roman"],R y1,{-1,0}],Text[Style["(a)",20,FontFamily->"Times New Roman"],{0.5L,-0.8L,L}]},Lighting->"Neutral",ViewPoint->{10,10,10},Boxed->False,PlotRange->L{{-1,1},{-1,1},{-1,1}}]


cy1=y1 . y y /Norm[y]^2


rady1=Norm[cy1-y1]


cx1=x1 . y y/Norm[y]^2


radx1=Norm[cx1-x1]


plb=Graphics3D[{{Directive[Opacity[0],EdgeForm[{Dashed,Black}]],ResourceFunction["Disk3D"][{0,0,0},1,y]},{Directive[White,Opacity[0.7],EdgeForm[{Blue}]],ResourceFunction["Disk3D"][cy1,rady1,y]},{Directive[Opacity[0],EdgeForm[{Red}]],ResourceFunction["Disk3D"][cx1,radx1,y]},{Red,Arrow[{{0,0,0},R x}],Arrow[{{0,0,0},R y}],Arrow[{{0,0,0},R z}]},{Blue,Arrow[{{0,0,0},R x1}],Arrow[{{0,0,0},R y1}],Arrow[{{0,0,0},R z1}]},{Darker[Green],Arrow[{{0,0,0},R x2}],Arrow[{{0,0,0},R y2}],Arrow[{{0,0,0},R z2}]},{Green,Arrowheads[0.03],Arrow[Tube[{x1,x2},0.015]],Arrow[Tube[{z1,z2},0.015]]},{Black,Arrowheads[0.03],Arrow[BezierCurve[Reverse[Table[0.8 R y+0.2 {Cos[\[Theta]],0,Sin[\[Theta]]},{\[Theta],-0.4,\[Pi],0.1}]]],-0.1]},Text[Style["\!\(\*SubscriptBox[\(y\), \(1\)]\)",20,FontFamily->"Times New Roman"],R y1,{-1,0}],Text[Style["\!\(\*SubscriptBox[\(y\), \(2\)]\)",20,FontFamily->"Times New Roman"],R y2,{-1,0}],Text[Style["\!\(\*SubscriptBox[\(x\), \(1\)]\)",20,FontFamily->"Times New Roman",Background->White],R x1,{1,0}],Text[Style["\!\(\*SubscriptBox[\(x\), \(2\)]\)",20,FontFamily->"Times New Roman"],R x2,{0,1}],Text[Style["y",20,FontFamily->"Times New Roman"],R y,{-1,0}],Text[Style["\!\(\*SubscriptBox[\(z\), \(1\)]\)=\!\(\*
StyleBox[\"z\",\nFontSlant->\"Italic\"]\)",20,FontFamily->"Times New Roman"],R z1,{0,-1}],Text[Style["\!\(\*SubscriptBox[\(z\), \(2\)]\)",20,FontFamily->"Times New Roman"],R z2,{0,-1}],Text[Style["(b)",20,FontFamily->"Times New Roman"],{0.5L,-0.8L,L}]},Lighting->"Neutral",ViewPoint->{10,10,6},Boxed->False,PlotRange->L{{-1,1},{-1,1},{-1,1}}]


cy2=y2 . z z 


rady2=Norm[cy2-y2]


cx2=x2 . z z


radx2=Norm[cx2-x2]


cz2=z2 . z z


radz2=Norm[cz2-z2]


plc=Graphics3D[{{Directive[Opacity[0],EdgeForm[{Dashed,Black}]],ResourceFunction["Disk3D"][cz2,radz2,z]},{Directive[White,Opacity[0.7],EdgeForm[{Blue}]],ResourceFunction["Disk3D"][cy2,rady2,z]},{Directive[Opacity[0],EdgeForm[{Red}]],ResourceFunction["Disk3D"][cx2,radx2,z]},{Red,Arrow[{{0,0,0},R x}],Arrow[{{0,0,0},R y}],Arrow[{{0,0,0},R z}]},(*{Blue,Arrow[{{0,0,0},R x1}],Arrow[{{0,0,0},R y1}],Arrow[{{0,0,0},R z1}]},*){Darker[Green],Arrow[{{0,0,0},R x2}],Arrow[{{0,0,0},R y2}],Arrow[{{0,0,0},R z2}]},{Purple,Arrow[{{0,0,0},R x3}],Arrow[{{0,0,0},R y3}],Arrow[{{0,0,0},R z3}]},{Purple,Arrowheads[0.03],Arrow[Tube[{x2,x3},0.015]],Arrow[Tube[{y2,y3},0.015]],Arrow[Tube[{z2,z3},0.01]]},{Black,Arrowheads[0.03],Arrow[BezierCurve[Table[0.8 R z+0.2 {Cos[\[Theta]],Sin[\[Theta]],0},{\[Theta],-0.4,\[Pi],0.1}]],-0.1]},(*Text[Style["Subscript[y, 1]",20,FontFamily->"Times New Roman"],R y1,{-1,0}],*)Text[Style["\!\(\*SubscriptBox[\(y\), \(2\)]\)",20,FontFamily->"Times New Roman"],R y2,{-1,0}],Text[Style["\!\(\*
StyleBox[\"y\",\nFontSlant->\"Italic\"]\)'",20,FontFamily->"Times New Roman"],R y3,{-1,0}],(*Text[Style["Subscript[x, 1]",20,FontFamily->"Times New Roman",Background->White],R x1,{1,0}],*)Text[Style["\!\(\*SubscriptBox[\(x\), \(2\)]\)",20,FontFamily->"Times New Roman"],R x2,{0,1}],Text[Style["\!\(\*
StyleBox[\"x\",\nFontSlant->\"Italic\"]\)'",20,FontFamily->"Times New Roman"],R x3,{0,1}],Text[Style["\!\(\*
StyleBox[\"z\",\nFontSlant->\"Italic\"]\)",20,FontFamily->"Times New Roman"],R z1,{0,-1}],Text[Style["\!\(\*SubscriptBox[\(z\), \(2\)]\)",20,FontFamily->"Times New Roman"],R z2,{1,-1}],Text[Style["\!\(\*
StyleBox[\"z\",\nFontSlant->\"Italic\"]\)'",20,FontFamily->"Times New Roman",Background->White],R z3,{-2.5,1}],Text[Style["(c)",20,FontFamily->"Times New Roman"],{0.5L,-0.8L,L}]},Lighting->"Neutral",ViewPoint->{10,10,6},Boxed->False,PlotRange->L{{-1,1},{-1,1},{-1,1}}]


Export["fig1_4.pdf",GraphicsRow[{pla,plb,plc},0,ImageSize->1024]]


{x1,y1,z1}=(r1=RotationTransform[65 degrad,z])/@{x,y,z};


{x2,y2,z2}=(r2=RotationTransform[25 degrad,y1])/@{x1,y1,z1};


{x3,y3,z3}=(r2=RotationTransform[25 degrad,z2])/@{x2,y2,z2};


Clear[arc3d];arc3d[st_,en_,n_:10]:=Block[{\[Delta]=en-Projection[en,st] },Print[\[Delta]];Table[Normalize[st (n-1-i)/(n-1)+en i/(n-1)] Norm[en] ,{i,0,n-1}]]


g15=Graphics3D[{{Directive[White,Opacity[0.7],EdgeForm[{Black}]],ResourceFunction["Disk3D"][{0,0,0},1,z]},{Directive[Lighter[Red],Opacity[0.05],EdgeForm[{Red}]],ResourceFunction["Disk3D"][{0,0,0},1,z3]},{Red,Arrow[{{0,0,0},R x}],Arrow[{{0,0,0},R y}],Arrow[{{0,0,0},R z}]},(*{Blue,Arrow[{{0,0,0},R x1}],Arrow[{{0,0,0},R y1}],Arrow[{{0,0,0},R z1}]},*)
{Dashed,Line[{-y1,y1}]},{Purple,Arrow[{{0,0,0},R x3}],Arrow[{{0,0,0},R y3}],Arrow[{{0,0,0},R z3}]},{Darker[Green],Arrowheads[0.03],Arrow[Tube[arc3d[0.6y,0.6 y1,4],0.01]],Arrow[arc3d[0.9y1,0.9 y3,4]],Arrow[arc3d[1.z,1. z3,4]]},(*Text[Style["Subscript[y, 1]",20,FontFamily->"Times New Roman"],R y1,{-1,0}],*)Text[Style["\!\(\*
StyleBox[\"y\",\nFontSlant->\"Italic\"]\)",20,FontFamily->"Times New Roman"],R y,{-1,0}],Text[Style["\!\(\*
StyleBox[\"y\",\nFontSlant->\"Italic\"]\)'",20,FontFamily->"Times New Roman"],R y3,{-1,0}],(*Text[Style["Subscript[x, 1]",20,FontFamily->"Times New Roman",Background->White],R x1,{1,0}],*)Text[Style["\!\(\*
StyleBox[\"x\",\nFontSlant->\"Italic\"]\)",20,FontFamily->"Times New Roman"],R x,{0,1}],Text[Style["\!\(\*
StyleBox[\"x\",\nFontSlant->\"Italic\"]\)'",20,FontFamily->"Times New Roman"],R x3,{0,1}],Text[Style["\!\(\*
StyleBox[\"z\",\nFontSlant->\"Italic\"]\)",20,FontFamily->"Times New Roman"],R z,{0,-1}],Text[Style["\!\(\*
StyleBox[\"z\",\nFontSlant->\"Italic\"]\)'",20,FontFamily->"Times New Roman",Background->White],R z3,{-2.5,1}],Text[Style["\[Alpha]",20,FontFamily->"Times New Roman"],(y+y1)/2 0.6,{-3,0}],Text[Style["\[Beta]",20,FontFamily->"Times New Roman"],(z+z3)/2 1.0,{3,0}],Text[Style["\[Gamma]",20,FontFamily->"Times New Roman",Background->White],(y1+y3)/2 0.9,{2,0}]},Lighting->"Neutral",ViewPoint->{10,10,6},Boxed->False,PlotRange->L{{-1,1},{-1,1},{-1,1}}]


Export["fig1_5.pdf",g15,ImageSize->1024]


nv=Normalize[{-0.7,0.3,0.7}];\[Theta]=50 degrad;{xp,yp,zp}=(RotationTransform[\[Theta],nv])/@{x,y,z};
xn=Projection[x,nv];yn=Projection[y,nv];zn=Projection[z,nv];
rx=Norm[xn-x];ry=Norm[yn-y];rz=Norm[zn-z];


\[CapitalTheta]=ArcTan[nv[[3]],Sqrt[nv[[1]]^2+nv[[2]]^2]]


\[CapitalPhi]=ArcTan[nv[[1]],nv[[2]]]


g16=Graphics3D[{Arrow[{-R nv,R nv}],{White,Directive[Opacity[0.5],EdgeForm[{Black}]],ResourceFunction["Disk3D"][xn,rx,nv]},{Black,Arrowheads[0.03],Arrow[BezierCurve[Table[0.8 R nv+0.2RotationTransform[\[Theta],nv][{1,0,0}],{\[Theta],-\[Pi]/2,3\[Pi]/4,0.1}]],-0.1]},{Directive[White,Opacity[0.7],EdgeForm[{Blue}]],ResourceFunction["Disk3D"][yn,ry,nv]},{Directive[Opacity[0],EdgeForm[{Red}]],ResourceFunction["Disk3D"][zn,rz,nv]},{Red,Arrow[{{0,0,0},R x}],Arrow[{{0,0,0},R y}],Arrow[{{0,0,0},R z}]},{Blue,Arrow[{{0,0,0},R xp}],Arrow[{{0,0,0},R yp}],Arrow[{{0,0,0},R zp}]},Text[Style["\!\(\*
StyleBox[\"y\",\nFontSlant->\"Italic\"]\)",20,FontFamily->"Times New Roman"],R y,{-1,0}],Text[Style["\!\(\*
StyleBox[\"y\",\nFontSlant->\"Italic\"]\)'",20,FontFamily->"Times New Roman"],R yp,{-1,0}],(*Text[Style["Subscript[x, 1]",20,FontFamily->"Times New Roman",Background->White],R x1,{1,0}],*)Text[Style["\!\(\*
StyleBox[\"x\",\nFontSlant->\"Italic\"]\)",20,FontFamily->"Times New Roman"],R x,{0,1}],Text[Style["\!\(\*
StyleBox[\"x\",\nFontSlant->\"Italic\"]\)'",20,FontFamily->"Times New Roman"],R xp,{0,1}],Text[Style["\!\(\*
StyleBox[\"z\",\nFontSlant->\"Italic\"]\)",20,FontFamily->"Times New Roman"],R z,{0,-1}],Text[Style["\!\(\*
StyleBox[\"z\",\nFontSlant->\"Italic\"]\)'",20,FontFamily->"Times New Roman"],R zp,{0,-1}],Text[Style["\!\(\*
StyleBox[\"n\",\nFontWeight->\"Plain\",\nFontSlant->\"Italic\"]\)",20,FontFamily->"Times New Roman"],R nv,{0,-1}],Text[Style["\[Omega]",20,FontFamily->"Times New Roman"],R 0.6nv,{-1.5,-0.5}]},Lighting->"Neutral",ViewPoint->{10,2,12},Boxed->False,PlotRange->L{{-1,1},{-1,1},{-1,1}}]


Export["fig1_6.pdf",g16,ImageSize->1024]


g17=Graphics3D[{{Red,Arrow[{{0,0,0},R x}],Arrow[{{0,0,0},R y}],Arrow[{{0,0,0},R z}]},{Green,Arrow[{{0,0,-1/2},R x+{0,0,-1/2}}],Arrow[{{0,0,-1/2},R y+{0,0,-1/2}}]},{Blue,Opacity[0.4],Sphere[{0,0,0},1/2]},{Directive[Opacity[0.0],EdgeForm[{Pink}]],ResourceFunction["Disk3D"][{0,0,0},1/2,{0,0,1}]},
{Pink,PointSize[Medium],Point[{1/3,1/3,1/6}]},
{Purple,Line[{{0,0,1/2},{1,1,-1/2}}]},
{Dashed,Line[{{1,0,-1/2},{1,1,-1/2},{0,1,-1/2}}]},Text[Style["\!\(\*
StyleBox[\"y\",\nFontSlant->\"Italic\"]\)",20,FontFamily->"Times New Roman"],R y,{-1,0}],Text[Style["\[Xi]",20,FontFamily->"Times New Roman"],R x+{0,0,-1/2},{0,1}],(*Text[Style["Subscript[x, 1]",20,FontFamily->"Times New Roman",Background->White],R x1,{1,0}],*)Text[Style["\!\(\*
StyleBox[\"x\",\nFontSlant->\"Italic\"]\)",20,FontFamily->"Times New Roman"],R x,{1,0}],Text[Style["\[Eta]",20,FontFamily->"Times New Roman"],R y+{0,0,-1/2},{0,1}],Text[Style["\!\(\*
StyleBox[\"z\",\nFontSlant->\"Italic\"]\)",20,FontFamily->"Times New Roman"],R z,{0,-1}],Text[Style["\[Zeta]",20,FontFamily->"Times New Roman"],{1,1,-1/2},{0,1}]},Lighting->"Neutral",ViewPoint->{10,2,5},Boxed->False,PlotRange->L{{-1/2,1},{-1/2,1},{-1,1}}]


Export["fig1_7.pdf",g17,ImageSize->1024]