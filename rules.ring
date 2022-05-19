define composite atom Si, Al, F, Cl, Re
define composite atom Zeo (heterogeneous site)
define composite atom HTA

input reactant "{Si}(O)"
input reactant "O({Al})({Si})"
input reactant "[{Al}+]O{Si}O{Si}"
input reactant "[{Re}](=O)(=O)(=O)(C)"
input reactant "[{F}H]"
input reactant "[{Cl}H]" 
input reactant "[H+]"
input reactant "[HH]"
input reactant "CC(=O)OC(=O)C" 				//acetic anhydride
input reactant "CC(=O)O" 					//acetic acid
input reactant "CC1=CC=CO1" 				//2-methylfuran
input reactant "CC1=CC=C(O1)C" 				//2,5-dimethylfuran 
input reactant "C1=CC=C(O1)" 				//furan
//input reactant "C(C1C(C(C(C(O1)O)O)O)O)O" 	//glucose 
//input reactant "O1CCCC1" 					//tetrahydrofuran  
//input reactant "CC1CCCO1" 					//2-methyltetrahydrofuran
//input reactant "CC1CCC(O1)C" 				//2,5-dimethyltetrahydrofuran 

global constraints on mol {
   fragment f1 {
      C labeled c1
	  C labeled c2 double bond to c1
	  C labeled c3 double bond to c1
   }
   fragment f2 {
      C labeled c1
	  C labeled c2 double bond to c1
	  O labeled o1 double bond to c2
   }
   fragment f3 {
      C labeled c1
	  C labeled c2 double bond to c1
	  O+ labeled o1 double bond to c2
   }
   ! mol contains f1 && ! mol contains f2 && ! mol contains f3
   mol.size <= 7
}

//acylation
rule aceticacid {
   reactant r1 { 
      C labeled c1
	  O labeled o1 double bond to c1
	  O labeled o2 single bond to c1
	  H labeled h1 single bond to o2
	  C labeled c2 single bond to c1
   }
   reactant r2 {
      O labeled o3
	  Si labeled s1 single bond to o3 
	  H labeled h2 single bond to o3
   }
   constraints {
      r2.size <= 2 && r1.size <= 4							
   }
   break bond (o3, h2)
   break bond (c1, o2)
   form bond (o2, h2)
   modify atomtype (c1, C+)
   modify atomtype (o3, O-)
}

rule aceticanhyd1 {
   reactant r1 {
      C labeled c1 {connected to C with single bond}
	  O labeled o1 double bond to c1
	  O labeled o2 single bond to c1
	  C labeled c2 single bond to o2 {connected to C with single bond}
	  O labeled o3 double bond to c2
   }
   reactant r2 {  
      Si labeled s1
	  O labeled o4 single bond to s1
	  H labeled h1 single bond to o4
   }
   constraints {
      fragment f1 {
	     O labeled o1 {in ring of size >= 4}
      }
	  fragment f2 {
	     C+ labeled c1
      }
	  ! r1 contains f1 && ! r1 contains f2
	  r2.size <= 2
   }
   decrease bond order (c1, o1)
   break bond (o4, h1)
   form bond (o1, h1)
   modify atomtype (o4, O-)
   modify atomtype (c1, C+)
}

rule aceticanhyd2 {
   reactant r1 {
      C labeled c1 {connected to C with single bond}
	  O labeled o1 double bond to c1
	  O labeled o2 single bond to c1
	  C+ labeled c2 single bond to o2 {connected to C with single bond}
	  O labeled o3 single bond to c2
   }
   reactant r2 {
      O labeled o4
	  Si labeled s1 single bond to o4
	  Al labeled a1 single bond to o4
   }
   constraints {
      r2.size <= 3
   }
   break bond (o2, c1)
   increase bond order (o2, c2)
   modify atomtype (c2, C)
   break bond (o4, a1)
   modify atomtype (a1, Al+)
   form bond (o4, c1)
}

rule SiOAl {
   reactant r1 {
      Al+ labeled a1
   }
   reactant r2 {
      Si labeled s1
	  O- labeled o1 single bond to s1
   }
   constraints {
      r1.size <= 1 && r2.size <= 2
   }
   form bond (o1, a1)
   modify atomtype (o1, O)
   modify atomtype (a1, Al)
}

rule aceticanhyd3 {
   reactant r1 {
      O labeled o1
	  Si labeled s1 single bond to o1
	  Al labeled a1 single bond to o1
   }
   reactant r2 { 
      O labeled o2
	  Si labeled s2 single bond to o2
	  C labeled c1 single bond to o2
	  O labeled o3 double bond to c1
	  C labeled c2 single bond to c1
   }
   constraints {
      fragment f1 {
	     C+ labeled c1
      }
	  fragment f2 {
	     C labeled c1
		 C labeled c2 double bond to c1
      }
	  ! r2 contains f1 && ! r2 contains f2
	  r1.size <= 3 && r2.size <= 5
   }
   break bond (o2, c1)
   break bond (o1, a1)
   form bond (a1, o2)
   modify atomtype (c1, C+)
   modify atomtype (o1, O-)
}

rule whelandintmf {
   reactant r1 {
      C+ labeled c1
	  O labeled o1 double bond to c1
   }
   reactant r2 {
      c labeled c2
	  c labeled c3 aromatic bond to c2
   }
   constraints {
      fragment f1 {
	     c labeled c1
		 C labeled c2 single bond to c1
		 O labeled o1 double bond to c2
		 C labeled c3 single bond to c2 {connected to 3 H with single bond}
      }
	  ! r2 contains f1
   }
   modify bond (c2, c3, single)
   form bond (c1, c2)
   modify atomtype (c1, C)
   modify atomtype (c3, C+)
}   

rule whelandint {
   reactant r1 {
      C+ labeled c1
	  O labeled o1 double bond to c1
   }
   reactant r2 {
      C labeled c2 {in ring of size >= 4}
	  C labeled c3 double bond to c2 {in ring of size >= 4}
   } 
   constraints {
      fragment f1 {
	     C labeled c1 {in ring of size >= 4}
		 C labeled c2 single bond to c1
		 O labeled o1 double bond to c2
		 C labeled c3 single bond to c2 {connected to 3 H with single bond}
      }
	  fragment f2 {
	     C+ labeled c1
      }
	  ! r2 contains f1 && ! r2 contains f2
   }
   decrease bond order (c2, c3)
   form bond (c1, c2)
   modify atomtype (c1, C)
   modify atomtype (c3, C+)
}

rule acylationend {
   reactant r1 {
      C labeled c1
	  C+ labeled c2 single bond to c1
	  C labeled c3 single bond to c1
	  O labeled o1 double bond to c3
	  H labeled h1 single bond to c1
   }
   reactant r2 {
      Si labeled s1
	  O- labeled o2 single bond to s1
   }
   constraints {
      fragment f1 {
	     c labeled c1
		 C labeled c2 single bond to c1
		 O labeled o1 double bond to c2
		 C labeled c3 single bond to c2 {connected to 3 H with single bond}
      }
	  fragment f2 {
	     C labeled c1 {in ring of size >= 4}
		 C labeled c2 single bond to c1
		 O labeled o1 double bond to c2
		 C labeled c3 single bond to c2 {connected to 3 H with single bond}
      }
	  r1 contains f1 || r1 contains f2 //only works on acylated molecules
   }
   break bond (c1, h1)
   increase bond order (c1, c2)
   form bond (o2, h1)
   modify atomtype (c2, C)
   modify atomtype (o2, O)
}

//alkylation
rule breakhf {
   reactant r1 {
      C labeled c1 {! connected to C with single bond, ! connected to C+}
	  C labeled c2 double bond to c1
	  C labeled c3 single bond to c2
   }
   reactant r2 {
      H labeled h1
	  F labeled f1 single bond to h1
   }
   constraints {
      fragment f1 {
	     C+ labeled c1
      }
	  fragment f2 {
	     O+ labeled o1
      }
	  fragment f3 {
	     O- labeled o1
      }
	  fragment f4 {
	     Si labeled s1
      }
	  fragment f5 {
	     O labeled o1
      }
	  ! r1 contains f1 && ! r1 contains f2 && ! r1 contains f3 && ! r1 contains f4 && ! r1 contains f5
   }
   decrease bond order (c1, c2)
   break bond (h1, f1)
   form bond (c1, h1)
   modify atomtype (c2, C+)
   modify atomtype (f1, F-)
}

rule alkylation {
   neutral reactant r1 { 
      c labeled c1
	  H labeled h1 single bond to c1
	  c labeled c2 aromatic bond to c1 {! connected to C with single bond, ! connected to 3 C with aromatic bond}
	  //c labeled c3 aromatic bond to c2 {! connected to C with single bond, ! connected to 3 C with aromatic bond}
	  //c labeled c4 aromatic bond to c3 {! connected to C with single bond, ! connected to 3 C with aromatic bond}
	  //c labeled c5 aromatic bond to c4 {! connected to C with single bond, ! connected to 3 C with aromatic bond}
	  //c labeled c6 aromatic bond to c5 {! connected to C with single bond, ! connected to 3 C with aromatic bond}
   }
   reactant r2 {
      C+ labeled c7
   }
   constraints {
	  fragment f2 {
	     c labeled c1
		 C labeled c2 single bond to c1
      }
	  fragment f3 {
	     O labeled o1
      }
	  ! r1 contains f2 //prevents multiple alkylations on same molecules
	  ! r2 contains f3 //only alkyl groups
   }
   modify bond (c1, c2, double)
   form bond (c1, c7)
   break bond (h1, c1)
   modify atomtype (c7, C)
   modify atomtype (h1, H+)
}

/*
rule alkylationnaphthalene {
   reactant r1 {
      c labeled c1
	  H labeled h1 single bond to c1
	  c labeled c2 aromatic bond to c1 {! connected to C with single bond}
	  c labeled c3 aromatic bond to c2 {! connected to C with single bond}
	  c labeled c4 aromatic bond to c3 {! connected to C with single bond}
	  c labeled c5 aromatic bond to c4 {! connected to C with single bond}
	  c labeled c6 aromatic bond to c5 {! connected to C with single bond}
	  c labeled c7 aromatic bond to c6 {! connected to C with single bond}
	  c labeled c8 aromatic bond to c7 {! connected to C with single bond}
      c labeled c9 aromatic bond to c8 {! connected to C with single bond}
	  c labeled c10 aromatic bond to c9 {! connected to C with single bond}
   }
   reactant r2 {
      C+ labeled c11
   }
   modify bond (c1, c2, double)
   form bond (c1, c11)
   break bond (h1, c1)
   modify atomtype (c11, C)
   modify atomtype (h1, H+)
}*/

rule catregen {
   reactant r1 {
      H+ labeled h1
   }
   reactant r2 {
      F- labeled f1
   }
   form bond (h1, f1)
   modify atomtype (f1, F)
   modify atomtype (h1, H)
}

//HAA
rule breakhCl {
   reactant r1 {
      H labeled h1
	  Cl labeled cl1 single bond to h1
   }
   reactant r2 {
      C labeled c1 {!connected to >= 2 O with double bond}
	  O labeled o1 double bond to c1
   }
   constraints {
      fragment f1 {
	     C+ labeled c1
      }
	  fragment f2 {
	     O+ labeled o1
      }
	  fragment f3 {
	     O- labeled o1
      }
	  fragment f4 {
	     Si labeled s1
      }
	  ! r2 contains f1 && ! r2 contains f2 && ! r2 contains f3 && ! r2 contains f4 
   }
   decrease bond order (c1, o1)
   break bond (cl1, h1)
   form bond (o1, h1)
   modify atomtype (cl1, Cl-)
   modify atomtype (c1, C+)
}

rule HAA1 {
   reactant r1 {
      C+ labeled c1
	  O labeled o1 single bond to c1
	  H labeled h1 single bond to o1
   }
   reactant r2 {
      c labeled c2 {connected to o with aromatic bond, ! connected to C with single bond, ! connected to C+ with single bond}
	  c labeled c3 aromatic bond to c2 {! connected to C with single bond}
   }
   modify bond (c2, c3, single)
   form bond (c2, c1)
   modify atomtype (c1, C)
   modify atomtype (c3, C+)
}

rule HAAdehydration {
   reactant r1{
      O labeled o1
      H labeled h1 single bond to o1
      C labeled c1 single bond to o1 {connected to =3 C with single bond} //methyl carbon
      C labeled c2 single bond to c1 {in ring of size >=3}
      H labeled h2 single bond to c2
      C+ labeled c3 single bond to c2
   }
   break bond (c1, o1)
   break bond (c2, h2)
   form bond (o1,h2)
   modify atomtype (c3, C)
   modify atomtype (c1, C+)
   increase bond order (c2, c3)
}

rule alk {
   reactant r1 {
      C+ labeled c1 
	  C labeled c2 single bond to c1
	  C labeled c3 single bond to c1
	  c labeled c4 single bond to c1
   }
   reactant r2 {
      c labeled c5 {connected to o with aromatic bond, ! connected to C with single bond, ! connected to C+ with single bond}
	  c labeled c6 aromatic bond to c5 {! connected to C with single bond}
   }
   modify bond (c5, c6, single)
   form bond (c5, c1)
   modify atomtype (c1, C)
   modify atomtype (c6, C+)
}   

rule HAAdeprotonate {
   reactant r1 {
      Cl- labeled cl1
   }
   reactant r2 {
      C labeled c1 {in ring of size >=4, connected to 2 C with single bond}
	  H labeled h1 single bond to c1
	  C+ labeled c2 single bond to c1 {in ring of size >= 4, connected to 2 C with single bond, connected to 2 H with single bond}
   }
   constraints {
      fragment f1 {
	     C labeled c1 {in ring of size >= 4}
		 C labeled c2 single bond to c1
		 O labeled o1 double bond to c2
		 C labeled c3 single bond to c2 {connected to 3 H with single bond}
      }
	  fragment f2 {
	     C labeled c1 {in ring of size >= 4}
		 C labeled c2 single bond to c1
		 O labeled o1 double bond to c2
		 O labeled o2 single bond to c2
      }
	  ! r2 contains f1 && ! r2 contains f2
   }
   break bond (c1, h1)
   form bond (cl1, h1)
   increase bond order (c1, c2)
   modify atomtype (cl1, Cl)
   modify atomtype (c2, C)
}

//ketonization
rule ketdehydration {
   reactant r1 {
      C labeled c1
	  O labeled o1 double bond to c1
	  O labeled o2 single bond to c1
	  H labeled h1 single bond to o2
   }
   reactant r2 {
      Si labeled s1
	  O labeled o3 single bond to s1
	  H labeled h2 single bond to o3
   }
   break bond (c1, o2)
   break bond (o3, h2)
   form bond (o2, h2)
   form bond (o3, c1)
}

rule zeoprotonation {
   reactant r1 {
      C labeled c1
	  O labeled o1 double bond to c1
	  O labeled o2 single bond to c1
	  H labeled h1 single bond to o2
   }
   reactant r2 {
      Al+ labeled a1
      O labeled o3 single bond to a1
	  Si labeled s1 single bond to o3
	  O labeled o4 single bond to s1
	  Si labeled s2 single bond to o4
   }
   constraints {
      fragment f1 {
	     C+ labeled c1
      }
	  fragment f2 {
	     O+ labeled o1
      }
	  fragment f3 {
	     O- labeled o1
      }
	  ! r1 contains f1 && ! r1 contains f2 && ! r1 contains f3
   }
   break bond (o4, s1)
   modify atomtype (s1, Si+)
   break bond (o2, h1)
   modify atomtype (o2, O-)
   form bond (o4, h1)
}

rule ketonization { 
   reactant r1 {
      Si labeled s1
	  O labeled o1 single bond to s1
	  C labeled c1 single bond to o1 {connected to O with double bond, connected to C with single bond}
   }
   reactant r2 {
      C labeled c2 {connected to O with double bond}
	  O- labeled o2 single bond to c2
	  C labeled c3 single bond to c2
   }
   constraints {
      fragment f1 {
	     C+ labeled c1
      }
	  fragment f2 {
	     O+ labeled o1
      }
	  fragment f3 {
	     O- labeled o1
      }
	  ! r1 contains f1 && ! r1 contains f2 && ! r1 contains f3
   }
   increase bond order (c2, o2)
   modify atomtype (o2, O)
   break bond (c2, c3)
   break bond (o1, c1)
   modify atomtype (o1, O-)
   form bond (c1, c3)
}

rule zeoregen {
   reactant r1 {
      Al+ labeled a1
	  O labeled o1 single bond to a1
	  Si+ labeled s1 single bond to o1
   }
   reactant r2 {
      Si labeled s2 
	  O- labeled o2 single bond to s2
   }
   constraints {
      r2.size <= 2
   }
   form bond (o2, a1)
   modify atomtype (a1, Al)
   modify atomtype (o2, O)
}

rule zeoregen2 {
   reactant r1 {
      Si labeled s1
	  O labeled o1 single bond to s1
	  Al labeled a1 single bond to o1
	  O labeled o2 single bond to a1
	  Si+ labeled s2 single bond to o2
   }
   reactant r2 {
      Si labeled s3
	  O labeled o3 single bond to s3
	  H labeled h1 single bond to o3
   }
   constraints {
      r1.size <= 5 && r2.size <= 2
   }
   break bond (o3, h1)
   break bond (o1, a1)
   form bond (o1, h1)
   form bond (s2, o3)
   modify atomtype (s2, Si)
   modify atomtype (a1, Al+)
}

//aldol
rule enol {
   reactant r1 {
      Cl labeled cl1
	  H labeled h1 single bond to cl1
   }
   reactant r2 {
      C labeled c1 {! connected to C with double bond}
	  O labeled o1 double bond to c1
	  C labeled c2 single bond to c1 
	  H labeled h2 single bond to c2
   }
   constraints {
      fragment f1 {
	     C+ labeled c1
      }
	  fragment f2 {
	     O+ labeled o1
      }
	  fragment f3 {
	     O- labeled o1
      }
	  fragment f4 {
	     Si labeled s1
      }
      fragment f5 {
	     O labeled o1
	     C labeled c1 double bond to o1
	     C labeled c2 single bond to c1
	     C labeled c3 single bond to c2
	     O labeled o2 single bond to c3
	     H labeled h1 single bond to o2
      }
	  fragment f6 {
         O labeled o1
	     C labeled c1 double bond to o1
	     C labeled c2 single bond to c1
	     C labeled c3 double bond to c2
      }
	  ! r2 contains f1 && ! r2 contains f2 && ! r2 contains f3 && ! r2 contains f4
	  ! r2 contains f5
	  ! r2 contains f6
   } //otherwise this rule takes the aldol product/dehydrated product as a reactant
   decrease bond order (c1, o1)
   break bond (h1, cl1)
   form bond (o1, h1)
   break bond (c2, h2)
   form bond (cl1, h2)
   increase bond order (c1, c2)
}

rule aldehydeprotonate {
   reactant r1 {
      Cl labeled cl1
	  H labeled h1 single bond to cl1
   }
   reactant r2 {
      C labeled c1 {! connected to O with single bond}
	  C labeled c2 single bond to c1
	  O labeled o1 double bond to c1
	  //H labeled h2 single bond to c1
   }
   constraints {
      fragment f1 {
	     C+ labeled c1
      }
	  fragment f2 {
	     O+ labeled o1
      }
	  fragment f3 {
	     O- labeled o1
      }
	  fragment f4 {
	     Si labeled s1
      }
      fragment f5 {
	     O labeled o1
	     C labeled c1 double bond to o1
	     C labeled c2 single bond to c1
	     C labeled c3 single bond to c2
	     O labeled o2 single bond to c3
	     H labeled h1 single bond to o2
      }
	  fragment f6 {
         O labeled o1
	     C labeled c1 double bond to o1
	     nonaromatic C labeled c2 single bond to c1
	     nonaromatic C labeled c3 double bond to c2
      }
	  ! r2 contains f1 && ! r2 contains f2 && ! r2 contains f3 && ! r2 contains f4
	  ! r2 contains f5
	  ! r2 contains f6
   }
   break bond (cl1, h1)
   form bond (o1, h1)
   modify atomtype (o1, O+)
   modify atomtype (cl1, Cl-)
}
	  
rule aldolcatregen {
   reactant r1 {
      Cl- labeled cl1
   }
   reactant r2 {
      nonaromatic C labeled c1 {! connected to >= 2 O with single bond}
	  C labeled c2 double bond to c1
	  O labeled o1 single bond to c1
	  H labeled h1 single bond to o1
   }
   constraints {
      fragment f1 {
	     C+ labeled c1
      }
	  fragment f2 {
	     O+ labeled o1
      }
	  fragment f3 {
	     O- labeled o1
      }
	  fragment f4 {
	     Si labeled s1
      }
	  ! r2 contains f1 && ! r2 contains f2 && ! r2 contains f3 && ! r2 contains f4 
   }
   break bond (o1, h1)
   form bond (cl1, h1)
   modify atomtype (cl1, Cl)
   modify atomtype (o1, O-)
}

rule aldol {
   reactant r1 {
      O- labeled o1
	  C labeled c1 single bond to o1
	  C labeled c2 double bond to c1
   }
   reactant r2 {
      O+ labeled o2
	  H labeled h1 single bond to o2
	  C labeled c3 double bond to o2
   }
   increase bond order (o1, c1)
   modify atomtype (o1, O)
   decrease bond order (c1, c2)
   decrease bond order (o2, c3)
   modify atomtype (o2, O)
   form bond (c2, c3)
}

rule aldoldehydration {
   reactant r1 {
      O labeled o1
	  C labeled c1 double bond to o1
	  C labeled c2 single bond to c1
	  H labeled h1 single bond to c2
	  C labeled c3 single bond to c2
	  O labeled o2 single bond to c3
	  H labeled h2 single bond to o2
   }
   reactant r2 {
      H labeled h3
	  Cl labeled cl1 single bond to h3
   }
   break bond (h3, cl1)
   form bond (o2, h3)
   break bond (c2, h1)
   form bond (cl1, h1)
   break bond (c3, o2)
   increase bond order (c2, c3)
}

//diels-alder
rule dielsalder {
   reactant r1 {
      C labeled c1
	  C labeled c2 double bond to c1
   }
   reactant r2 {
      C labeled c3
	  C labeled c4 double bond to c3
	  C labeled c5 single bond to c4
	  C labeled c6 double bond to c5
   }
   constraints {
      fragment f {
         C labeled c3
	     C labeled c4 double bond to c3
	     C labeled c5 single bond to c4
	     C labeled c6 double bond to c5
      }
      fragment f1 {
	     C+ labeled c1
      }
	  fragment f2 {
	     O+ labeled o1
      }
	  fragment f3 {
	     O- labeled o1
      }
	  fragment f4 {
	     Si labeled s1
      }
	  ! r1 contains f
	  ! r1 contains f1 && ! r1 contains f2 && ! r1 contains f3 && ! r1 contains f4
	  ! r2 contains f1 && ! r2 contains f2 && ! r2 contains f3 && ! r2 contains f4
   }
   decrease bond order (c1, c2)
   decrease bond order (c3, c4)
   decrease bond order (c5, c6)
   increase bond order (c4, c5)
   form bond (c1, c3)
   form bond (c2, c6)
}  

// for 2,5-dimethylfuran
rule dielsalder2 {
   reactant r1 {
      C labeled c1
	  C labeled c2 double bond to c1
   }
   reactant r2 {
	  c labeled c3 //{! connected to 2 c with aromatic bond}
	  c labeled c4 aromatic bond to c3
	  c labeled c5 aromatic bond to c4
	  c labeled c6 aromatic bond to c5
   }
   constraints {
      fragment f {
         C labeled c3
	     C labeled c4 double bond to c3
	     C labeled c5 single bond to c4
	     C labeled c6 double bond to c5
      }
	  ! r1 contains f
   }
   decrease bond order (c1, c2)
   modify bond (c3, c4, single)
   modify bond (c4, c5, double)
   modify bond (c5, c6, single)
   form bond (c1, c3)
   form bond (c2, c6)
}

rule dadehydration {
   reactant r1 {
      C labeled c1 {connected to >= 2 C with single bond}
	  C labeled c2 single bond to c1 {! connected to C with double bond, in ring of size > 4}
	  H labeled h1 single bond to c2
	  O labeled o1 single bond to c1
	  C labeled c3 single bond to o1 {connected to >= 2 C with single bond}
	  C labeled c4 single bond to c3 {! connected to C with double bond, in ring of size > 4}
	  H labeled h2 single bond to c4
   }
   constraints {
      fragment f1 {
	     C+ labeled c1
      }
	  fragment f2 {
	     O+ labeled o1
      }
	  fragment f3 {
	     O- labeled o1
      }
	  fragment f4 {
	     Si labeled s1
	  }
	  ! r1 contains f1 && ! r1 contains f2 && ! r1 contains f3 && ! r1 contains f4
   }
   break bond (c2, h1)
   break bond (c1, o1)
   form bond (o1, h1)
   increase bond order (c1, c2)
   break bond (c4, h2)
   break bond (c3, o1)
   form bond (o1, h2)
   increase bond order (c3, c4)
}  

//polyol dehydration
rule alcdehyd1 {
   reactant r1 {
      C labeled c1
	  O labeled o1 single bond to c1
	  H labeled h1 single bond to o1
	  C labeled c2 single bond to c1 {! connected to O with single bond}
	  C labeled c3 single bond to c1 {! connected to O with single bond}
	  H labeled h2 single bond to c3
   }
   reactant r2 {
      Re labeled re1 {connected to 3 O with double bond}
	  O labeled o2 double bond to re1
	  O labeled o3 double bond to re1
   }
   constraints {
      fragment f1 {
         C+ labeled c1
      }
	  fragment f2 {
	     O- labeled o1
      }
	  fragment f3 {
	     O+ labeled o1
      }
	  fragment f4 {
	     Si labeled s1
      }
	  fragment f5 {
	     O labeled o1
	     C labeled c1 double bond to o1
	     C labeled c2 single bond to c1
	     C labeled c3 single bond to c2
	     O labeled o2 single bond to c3
	     H labeled h1 single bond to o2
      }
	  ! r1 contains f1 && ! r1 contains f2 && ! r1 contains f3 && ! r1 contains f4 && ! r1 contains f5
   }
   break bond (o1, h1)
   decrease bond order (re1, o2)
   form bond (o2, h1)
   form bond (o1, re1)
   break bond (c1, o1)
   modify atomtype (c1, C+)
   modify atomtype (o1, O-)
}

rule alcdehyd2 {
   reactant r1 {
      C+ labeled c1 {connected to H with single bond}
	  C labeled c2 single bond to c1
	  C labeled c3 single bond to c1
	  H labeled h1 single bond to c3
   }
   reactant r2 {
      Re labeled re1 {connected to 2 O with double bond, ! connected to >= 2 O with single bond, ! connected to >= 2 O-} 
	  O- labeled o1 single bond to re1
	  O labeled o2 single bond to re1
	  H labeled h2 single bond to o2
   }
   constraints {
      fragment f1 {
	     O labeled o1
		 H labeled h1 single bond to o1
      }
	  fragment f2 {
	     O+ labeled o1
      }
	  fragment f3 {
	     O- labeled o1
      }
	  fragment f4 {
	     Si labeled s1
	  }
	  fragment f5 {
	     O labeled o1
	     C labeled c1 double bond to o1
	     C labeled c2 single bond to c1
	     C labeled c3 single bond to c2
	     O labeled o2 single bond to c3
	     H labeled h1 single bond to o2
      }
	  ! r1 contains f1 && ! r1 contains f2 && ! r1 contains f3 && ! r1 contains f4 && ! r1 contains f5
   }
   break bond (c3, h1)
   break bond (re1, o2)
   form bond (o2, h1)
   increase bond order (re1, o1)
   modify atomtype (o1, O)
   increase bond order (c1, c3)
   modify atomtype (c1, C)
}

rule polyoldehyd1 {
   reactant r1 { 
      C labeled c1
	  O labeled o1 single bond to c1
	  H labeled h1 single bond to o1
	  C labeled c2 single bond to c1
	  O labeled o2 single bond to c2
	  H labeled h2 single bond to o2
   }
   reactant r2 {
      Re labeled re1 {connected to 3 O with double bond}
	  O labeled o3 double bond to re1
   }
   constraints {
      fragment f1 {
	     Re- labeled re1
      }
	  fragment f2 {
	     Re labeled re1
      }
	  ! r1 contains f1
	  ! r1 contains f2
   }
   break bond (o1, h1)
   decrease bond order (re1, o3)
   form bond (o3, h1)
   form bond (o1, re1)
   form bond (o2, re1)
   modify atomtype (o2, O+)
   modify atomtype (re1, Re-)
}

rule polyoldehyd2 {
   reactant r1 {
      Re- labeled re1
      O labeled o1 single bond to re1
	  H labeled h1 single bond to o1
	  O+ labeled o2 single bond to re1
	  H labeled h2 single bond to o2
	  O labeled o3 double bond to re1
   }
   reactant r2 {
      C labeled c1
	  O labeled o4 single bond to c1
	  H labeled h3 single bond to o4
	  C labeled c2 single bond to c1 {! connected to O with any bond, ! connected to C with double bond}
	  C labeled c3 single bond to c1 //might need to change this later
   }
   constraints {
      fragment f1 {
	     Re- labeled re1
      }
	  fragment f2 {
	     Re labeled re1
      }
	  ! r2 contains f1
	  ! r2 contains f2
   }
   break bond (o2, h2)
   modify atomtype (o2, O)
   break bond (re1, o1)
   modify atomtype (re1, Re)
   form bond (o1, h2)
   decrease bond order (re1, o3)
   break bond (o4, h3)
   form bond (o3, h3)
   form bond (o4, re1)
}

rule polyoldehyd3 {
   reactant r1 {
      Re labeled re1
	  O labeled o1 single bond to re1
	  H labeled h1 single bond to o1
	  O labeled o2 single bond to re1
	  C labeled c1 single bond to o2
	  H labeled h2 single bond to c1
	  C labeled c2 single bond to c1 {! connected to O with any bond}
	  O labeled o3 single bond to re1
	  C labeled c3 single bond to o3
	  C labeled c4 single bond to c3
	  O labeled o4 single bond to c4
	  ringbond o4 single bond to re1
   }
   break bond (re1, o1)
   break bond (c1, h2)
   form bond (o1, h2)
   break bond (re1, o2)
   increase bond order (c1, o2)
   break bond (o4, c4)
   break bond (o3, c3)
   increase bond order (re1, o3)
   increase bond order (re1, o4)
   increase bond order (c3, c4)
} 

//bronsted alcohol dehydration
rule bronstalk {
   reactant r1 {
      C labeled c1
	  O labeled o1 single bond to c1
	  H labeled h1 single bond to o1
	  C labeled c2 single bond to c1 {! connected to O with any bond}
   }
   reactant r2 {
      Si labeled s1
	  O labeled o2 single bond to s1
	  H labeled h2 single bond to o2
   }
   constraints {
      fragment f1 {
	     C+ labeled c1
      }
	  fragment f2 {
	     O+ labeled o1
      }
	  fragment f3 {
	     O- labeled o1
      }
	  fragment f4 {
	     Si labeled s1
	  }
	  fragment f5 {
	     O labeled o1
	     C labeled c1 double bond to o1
	     C labeled c2 single bond to c1
	     C labeled c3 single bond to c2
	     O labeled o2 single bond to c3
	     H labeled h1 single bond to o2
      }
	  ! r1 contains f1 && ! r1 contains f2 && ! r1 contains f3 && ! r1 contains f4 && ! r1 contains f5
   }
   break bond (o2, h2)
   break bond (c1, o1)
   form bond (o1, h2)
   form bond (o2, c1)
}

rule bronstalk2 {
   reactant r1 {
      Si labeled s1
	  O labeled o1 single bond to s1
	  C labeled c1 single bond to o1
	  C labeled c2 single bond to c1
	  H labeled h1 single bond to c2
   }
   constraints {
      fragment f1 {
	     C+ labeled c1
      }
	  fragment f2 {
	     O+ labeled o1
      }
	  fragment f3 {
	     O- labeled o1
      }
	  ! r1 contains f1 && ! r1 contains f2 && ! r1 contains f3 
   }
   break bond (c2, h1)
   break bond (c1, o1)
   increase bond order (c1, c2)
   form bond (o1, h1)
}

rule bronstether {
   reactant r1 {
      C labeled c1
	  O labeled o1 single bond to c1
	  H labeled h1 single bond to o1
	  C labeled c2 single bond to c1 {! connected to O with any bond}
   }
   reactant r2 {
      Si labeled s1
	  O labeled o2 single bond to s1
	  H labeled h2 single bond to o2
   }
   constraints {
      fragment f1 {
	     C+ labeled c1
      }
	  fragment f2 {
	     O+ labeled o1
      }
	  fragment f3 {
	     O- labeled o1
      }
	  fragment f4 {
	     Si labeled s1
	  }
	  fragment f5 {
	     O labeled o1
	     C labeled c1 double bond to o1
	     C labeled c2 single bond to c1
	     C labeled c3 single bond to c2
	     O labeled o2 single bond to c3
	     H labeled h1 single bond to o2
      }
	  ! r1 contains f1 && ! r1 contains f2 && ! r1 contains f3 && ! r1 contains f4 && ! r1 contains f5
   }
   break bond (o2, h2)
   break bond (c1, o1)
   form bond (o1, h2)
   modify atomtype (o2, O-)
   modify atomtype (c1, C+)
}

rule bronstether2 {
   reactant r1 {
      C labeled c1
	  O labeled o1 single bond to c1
	  H labeled h1 single bond to o1
	  C labeled c2 single bond to c1 {! connected to O with any bond}
   }
   reactant r2 {
      Si labeled s1
	  O- labeled o2 single bond to s1
   }
   constraints {
      fragment f1 {
	     C+ labeled c1
      }
	  fragment f2 {
	     O+ labeled o1
      }
	  fragment f3 {
	     O- labeled o1
      }
	  fragment f4 {
	     Si labeled s1
      }
	  fragment f5 {
	     O labeled o1
	     C labeled c1 double bond to o1
	     C labeled c2 single bond to c1
	     C labeled c3 single bond to c2
	     O labeled o2 single bond to c3
	     H labeled h1 single bond to o2
      }
	  ! r1 contains f1 && ! r1 contains f2 && ! r1 contains f3 && ! r1 contains f4 && ! r1 contains f5
   }
   break bond (o1, h1)
   form bond (o2, h1)
   modify atomtype (o2, O)
   modify atomtype (o1, O-)
}

rule bronstether3 {
   reactant r1 {
      C+ labeled c1
	  C labeled c2 single bond to c1
   }
   reactant r2 {
      O- labeled o1
	  C labeled c3 single bond to o1
   }
   constraints {
      fragment f1 {
	     C+ labeled c1
      }
	  fragment f2 {
	     O+ labeled o1
      }
	  fragment f3 {
	     O- labeled o1
      }
	  fragment f4 {
	     Si labeled s1
      }
	  ! r2 contains f1
	  ! r1 contains f2 && ! r2 contains f2
	  ! r1 contains f3
	  ! r1 contains f4 && ! r2 contains f4
   }
   form bond (o1, c1)
   modify atomtype (c1, C)
   modify atomtype (o1, O)
} 

//sugar isomerization
rule openglucose {
   reactant r1 {
      O labeled o1 {connected to 2 C, in ring of size >= 4}
	  C labeled c1 single bond to o1 {! connected to >=2 C, in ring of size >= 4}
	  O labeled o2 single bond to c1
	  H labeled h1 single bond to o2
   }
   reactant r2 {
      H+ labeled h2
   } 
   constraints {
      fragment f1 {
	     O labeled o1 {in ring of size >= 4}
		 C labeled c1 single bond to o1 {in ring of size >= 4}
		 O labeled o2 single bond to c1
		 C labeled c2 single bond to c1 {in ring of size >= 4}
		 O labeled o3 single bond to c2
		 C labeled c3 single bond to c2 {in ring of size >= 4}
		 O labeled o4 single bond to c3
		 C labeled c4 single bond to c3 {in ring of size >= 4}
		 O labeled o5 single bond to c4
      }
	  r1 contains f1
   }
   form bond (o1, h2)
   modify atomtype (h2, H)
   break bond (c1, o1)
   break bond (o2, h1)
   modify atomtype (h1, H+)
   increase bond order (c1, o2)
}  

rule sugarisomerization {
   reactant r1 {
      C labeled c1 {! connected to >= 2 C with single bond}
	  O labeled o1 double bond to c1
	  C labeled c2 single bond to c1
	  H labeled h1 single bond to c2
	  O labeled o2 single bond to c2
	  H labeled h2 single bond to o2
   }
   break bond (o2, h2)
   break bond (c2, h1)
   decrease bond order (c1, o1)
   form bond (c1, h1)
   increase bond order (c2, o2)
   form bond (o1, h2)
} 

rule closefructose {
   reactant r1 {
      C labeled c1 {connected to 2 C with single bond}
	  O labeled o1 double bond to c1
	  C labeled c2 single bond to c1
	  C labeled c3 single bond to c2
	  C labeled c4 single bond to c3
	  O labeled o2 single bond to c4
	  H labeled h1 single bond to o2
   }
   decrease bond order (c1, o1)
   form bond (c1, o2)
   break bond (o2, h1)
   form bond (o1, h1)
}

//hydrodeoxygenation
rule hdo1 { 
   reactant r1 {
      O labeled o1
	  C labeled c1 double bond to o1
   }
   reactant r2 {
      H labeled h1
	  H labeled h2 single bond to h1
   }
   constraints {
	   fragment f1 {
	      C+ labeled c1
	   }
	   fragment f2 {
	      O+ labeled o1
	   }
	   fragment f3 {
	      O- labeled o1
	   }
	   fragment f4 {
	      Si labeled s1
	   }
	   ! r1 contains f1 && ! r1 contains f2 && ! r1 contains f3 && ! r1 contains f4
   }
   decrease bond order (c1, o1)
   break bond (h1, h2)
   form bond (o1, h1)
   form bond (c1, h2)
}

rule hdo2 {
   reactant r1 {
      C labeled c1 
	  O labeled o1 single bond to c1
	  H labeled h1 single bond to o1
   }
   reactant r2 {
      H labeled h2
	  H labeled h3 single bond to h2
   }
   constraints {
      fragment f1 {
	     C+ labeled c1
      }
	  fragment f2 {
	     O+ labeled o1
      }
	  fragment f3 {
	     O- labeled o1
      }
	  fragment f4 {
	     Si labeled s1
	  }
	  fragment f5 {
	     C labeled c1
		 O labeled o1 double bond to c1
      }
	  ! r1 contains f1 && ! r1 contains f2 && ! r1 contains f3 && ! r1 contains f4 && ! r1 contains f5
   }
   form bond (o1, h2)
   break bond (h2, h3)
   break bond (c1, o1)
   form bond (c1, h3)
}

//dehydra-decyclization
rule ddc1 {
   reactant r1 {
      Si labeled s1
	  O labeled o1 single bond to s1
	  H labeled h1 single bond to o1
   }
   reactant r2 {
      O labeled o2 {in ring of size >=4}
	  C labeled c1 single bond to o2 
	  C labeled c2 single bond to c1 
	  H labeled h2 single bond to c2
	  C labeled c3 single bond to c2 
   }
   constraints {
      fragment f1 {
	     C+ labeled c1
      }
	  ! r2 contains f1
   }
   break bond (o1, h1)
   form bond (o2, h1)
   break bond (c1, o2)
   form bond (o1, c1)
}

rule ddc2 {
   reactant r1 {
      Si labeled s1
	  O labeled o1 single bond to s1
	  C labeled c1 single bond to o1 {! connected to > 1 C}
	  C labeled c2 single bond to c1
	  H labeled h1 single bond to c2
   }
   constraints {
      fragment f1 {
	     C+ labeled c1
      }
	  fragment f2 {
	     O+ labeled o1
      }
	  fragment f3 {
	     O- labeled o1
      }
	  fragment f4 {
	     O labeled o1
		 H labeled h1 single bond to o1
      }
	  ! r1 contains f1 && ! r1 contains f2 && ! r1 contains f3
	  r1 contains f4
   }
   break bond (c2, h1)
   form bond (c1, h1)
   break bond (c1, o1)
   form bond (o1, c2)
}

rule ddc3 {
   reactant r1 {
      Si labeled s1
	  O labeled o1 single bond to s1
	  C labeled c1 single bond to o1 {connected to >= 2 C}
	  C labeled c2 single bond to c1 {connected to >= 2 C}
	  H labeled h1 single bond to c2
   }
   constraints {
      fragment f1 {
	     C+ labeled c1
      }
	  fragment f2 {
	     O+ labeled o1
      }
	  fragment f3 {
	     O- labeled o1
      }
	  fragment f4 {
	     O labeled o1
		 H labeled h1 single bond to o1
      }
	  ! r1 contains f1 && ! r1 contains f2 && ! r1 contains f3
   }
   break bond (o1, c1)
   break bond (c2, h1)
   form bond (o1, h1)
   increase bond order (c1, c2)
}

rule ddc4 {
   reactant r1 {
      C labeled c1
	  H labeled h1 single bond to c1
	  C labeled c2 single bond to c1
	  C labeled c3 double bond to c2
	  C labeled c4 single bond to c3
	  O labeled o1 single bond to c4
   }
   constraints {
      fragment f1 {
	     C+ labeled c1
      }
	  fragment f2 {
	     O+ labeled o1
      }
	  fragment f3 {
	     O- labeled o1
      }
	  fragment f4 {
	     Si labeled s1
      }
	  fragment f5 {
	     O labeled o1
		 H labeled h1 single bond to o1
      }
	  ! r1 contains f1 && ! r1 contains f2 && ! r1 contains f3 && ! r1 contains f4
	  r1 contains f5
   }
   break bond (c1, h1)
   increase bond order (c1, c2)
   decrease bond order (c2, c3)
   increase bond order (c3, c4)
   break bond (c4, o1)
   form bond (o1, h1)
}

//hydrogenation
rule DoubleBondHydrogenation {
   reactant r1 {
      C labeled c1
	  C labeled c2 double bond to c1
   }
   reactant r2 {
      H labeled h1
	  H labeled h2 single bond to h1
   }
   constraints {
      fragment f1 {
	     C+ labeled c1
      }
	  fragment f2 {
	     O+ labeled o1
      }
	  fragment f3 {
	     O- labeled o1
      }
	  fragment f4 {
	     Si labeled s1
	  }
	  ! r1 contains f1 && ! r1 contains f2 && ! r1 contains f3 && ! r1 contains f4
   }
   break bond (h1, h2)
   decrease bond order (c1, c2)
   form bond (c1, h1)
   form bond (c2, h2)
}
rule AromaticHydrogenation {
   reactant r1 {
      aromatic $ labeled a1
	  aromatic $ labeled a2 aromatic bond to a1
   }
   reactant r2 {
      H labeled h1
	  H labeled h2 single bond to h1
   }
   modify bond (a1, a2, single)
   break bond (h1, h2)
   form bond (a1, h1)
   form bond (a2, h2)
}
rule FuranHydro {
   reactant r1{
      o labeled o1
      c labeled c1 aromatic bond to o1
      c labeled c2 aromatic bond to c1
   }
   reactant r2{
      H labeled h1
      H labeled h2 single bond to h1
   }
   break bond (h1,h2)
   modify bond (c1,c2,single)
   form bond (c1,h1)
   form bond (c2,h2)
} 

find all mol {
   fragment f1 {
      C+ labeled c1
   }
   fragment f2 {
      O+ labeled o1
   }
   fragment f3 {
      O- labeled o1
   }
   fragment f4 {
      Si labeled s1
   }
   fragment f5 {
      Re labeled r1
   }
   fragment f6 {
      C labeled c1
   }
   ! mol contains f1 && ! mol contains f2 && ! mol contains f3 && ! mol contains f4 && ! mol contains f5 
   mol contains f6
} store in "all-products.txt"

find pathways to mol {
   fragment f1 {
      C+ labeled c1
   }
   fragment f2 {
      O+ labeled o1
   }
   fragment f3 {
      O- labeled o1
   }
   fragment f4 {
      Si labeled s1
   }
   fragment f5 {
      Re labeled r1
   }
   fragment f6 {
      C labeled c1
   }
   ! mol contains f1 && ! mol contains f2 && ! mol contains f3 && ! mol contains f4 && ! mol contains f5 
   mol contains f6
} constraints {
   maximum length shortest + 0
   eliminate similar pathways
   contains 1 mol {mol is "o1cccc1"}
} store in "all-pathways.txt" 