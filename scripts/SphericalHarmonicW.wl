(* ::Package:: *)

(* ::Input:: *)
(*SphericalHarmonicW[J_,M_,\[CurlyTheta]_]:=1/(2(J+1)) ((J+M+1)SphericalHarmonicY[J+1/2,M+1/2,\[CurlyTheta],0]^2+(J-M+1)SphericalHarmonicY[J+1/2,M-1/2,\[CurlyTheta],0]^2)*)


(* ::Text:: *)
(*Test a relation from VKM's table:*)


(* ::Input:: *)
(*FullSimplify[SphericalHarmonicW[7/2,1/2,\[CurlyTheta]]-1/(64\[Pi]) (175 Cos[\[CurlyTheta]]^6-165 Cos[\[CurlyTheta]]^4+45 Cos[\[CurlyTheta]]^2+9)]*)


(* ::Input:: *)
(*SphericalHarmonicWperp[J_,M_,\[CurlyTheta]_]:=1/(2J(J+1)) ((J-M)(J+M+1)SphericalHarmonicY[J,M+1,\[CurlyTheta],0]^2+2M^2 SphericalHarmonicY[J,M+1,\[CurlyTheta],0]^2+(J+M)(J-M+1)SphericalHarmonicY[J,M-1,\[CurlyTheta],0]^2)*)


(* ::Input:: *)
(*SphericalHarmonicWpar[J_,M_,\[CurlyTheta]_]:=SphericalHarmonicY[J,M,\[CurlyTheta],0]^2*)


(* ::ChatDelimiter:: *)
(**)
