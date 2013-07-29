#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <vector>
#include <iomanip>
using namespace std;

double StoF(string CONV); 
int addatom_xyz(int CNT1, string LINE1); int addatom_mol(int CNT1, string LINE1); int addatom_gro(int CNT1, string LINE1); int addatom_pdb(int CNT1, string LINE1);
bool is_number(char c_a); int get_charge(string mol_sym);

string helpS="program made to convert geometry files into a dalton molecule input file. \
Requires 1 input files, a geometry in .mol, .gro, .xyz, or .pdb format. \
Usage: < -i {in.xyz} input> <-o {out.mol} > <-b basis set {6-31G*} > <-f file type(xyz,mol,gro,pdb)> <-n max read> \
<-q charge ><-c title/comment line> <-h print help and quit>\r\n\r\n\
for PDB format many different atom types may be possible. consider setting the '-p2' flag to only use the first character. \
This will make it impossible to use any atoms with more than 1 character symbols. Consider reformatting pdb if this is an issue.\r\n";
string COMMENT1="TITLE";
string mod_date="5-29-2013";
bool p2_flag=false;

struct ATOM_INF {
        int xy;
        //pubilc:
        float posx;
        float posy;
        float posz;
        string type;
} atom[12000];

int main (int argc, char *argv[])
{	 
cout << "Program 2gauss written by Cody Covington. Last Modified "<<mod_date<<". use -h for help. "<<endl;
if (argc < 2)
    {// If the user didn't provide a filename command line argument, print an error and exit.
        cout << "Usage: " << argv[0] << " < -i {in.xyz} input> <-o {out.mol} > <-d file to insert{gaussin.dat}> <-f file type(xyz,mol,gro,pdb)> <-n max read> " << endl<<" <-k checkpoint name> <-c title/comment line> <-h print help and quit> <-p2>"<<endl;
        return 1;
    }                 
string INFILE="in.xyz",OUTFILE="out.mol",basis="6-31G*",FTYPE=""; //default file in and out and data file
int MAXREAD=12000,i=1,NUM=0,a=0; 
double q_charge=0.0;
while (i + 1 <= argc){ // Check that we haven't finished parsing already
string ARG=argv[i];
if (ARG.find("-i")<1) {	INFILE=argv[i+1]; } 
else if (ARG.find("-o")<1) { OUTFILE=argv[i+1]; } 
else if (ARG.find("-b")<1) { basis=argv[i+1];  }
else if (ARG.find("-n")<1) { MAXREAD= (int)StoF(argv[i+1]); }
else if (ARG.find("-f")<1) { FTYPE=argv[i+1];   }
else if (ARG.find("-h")<1) { cout << helpS << endl; return 0; }
else if (ARG.find("-c")<1) { COMMENT1=argv[i+1];  }
else if (ARG.find("-p2")<1) { p2_flag = true; }
else if (ARG.find("-q")<1) { q_charge= StoF(argv[i+1]); }
i++;
}
string STRING;
ifstream infile;
    infile.open (INFILE.c_str());
if(infile.fail()==true){ cerr << "Input File "<<INFILE<<" Not Found" << endl; return 1; }
if(FTYPE=="mol"){ cout <<" input file .mol selected "<<endl; getline(infile,STRING); getline(infile,STRING); getline(infile,STRING); getline(infile,STRING);}
else if(FTYPE=="xyz"){ cout <<" input file .xyz selected "<<endl; getline(infile,STRING); getline(infile,STRING); } //read comment lines
else if(FTYPE=="gro"){ cout <<" input file .gro selected "<<endl; getline(infile,STRING); getline(infile,STRING); }
else if(FTYPE=="pdb"){
cout <<" input file PDB selected "<<endl;
//go until ATOM lines handeled by the addatom_pdb funciton
}
else { 
cerr <<" no Input file type given "<< "Usage: " << argv[0] << " < -i {in.mol} input> <-o {run.gjf} > <-d file to insert{gaussin.dat}> <-f file type(xyz,mol,gro)> <-n max read> " << endl; 
return 1; 
}
while(a<1) // To get you all the lines.
    { 
if(!infile.eof()==0){ a=1; cerr << "possible error in coordinate file" << endl; }
NUM++; 
   getline(infile,STRING); // Saves the line in STRING. //cout<<STRING<<endl;  
if( FTYPE=="mol" ){
   if(addatom_mol(NUM,STRING)==1){ NUM=NUM-1; a=1;} //end or error //else cout<<"line "<<NUM<<" read"<<endl;
}
if( FTYPE=="xyz" ){
   if(addatom_xyz(NUM,STRING)==1){ NUM=NUM-1; a=1;} //end or error //else cout<<"line "<<NUM<<" read"<<endl;
}
if( FTYPE=="gro" ){
   if(addatom_gro(NUM,STRING)==1){ NUM=NUM-1; a=1;} //end or error //else cout<<"line "<<NUM<<" read"<<endl;
}
if( FTYPE=="pdb" ){
   int b3=addatom_pdb(NUM,STRING);
   if(b3==1){ NUM=NUM-1; a=1;} //end or error //else cout<<"line "<<NUM<<" read"<<endl;
else if(b3==2){ 
NUM=NUM-1; //a row without "ATOM"
}
}

if(NUM>=MAXREAD) break;
    }
cout << NUM << " atoms read from "<< INFILE <<endl;

//print dalton crd file (.mol file)
//determine # of atomtypes. 
string atom_type[51]; //store different atom types
int atom_cnt[51] = {0}; int num_type=0;  // count number of atoms of each type and total number of different types
for(int jj=1;jj<=NUM;jj++){
 
for(int kk=1;kk<=50; kk++){
if(atom_type[kk] == ""){ //add new type
atom_type[kk] = atom[jj].type;
atom_cnt[kk]++; num_type++;
//cout <<" new atom type found "<<atom_type[kk]<< endl;
break;
}
if(atom[jj].type == atom_type[kk]){ //old type found
atom_cnt[kk]++;
//cout <<" atom type "<<atom_type[kk]<< " found atom #"<<atom_cnt[kk]<< endl;
break;
}
if(kk==50){
cerr<<"cannot handle coordinate file with more than 50 different atom types."<<endl;
return 1;
}
}
}
cout <<num_type<<" total different atom types found"<<endl;
//set up file
ofstream outfile2 ( OUTFILE.c_str() );
outfile2 << "ATOMBASIS"<<endl;
outfile2 <<"--COMMENT: "<< COMMENT1 <<endl;
outfile2 <<"--COMMENT: file created by 2dalton from crd file "<<INFILE<<endl;
//info line  Atomtypes=2 Charge=0.0 Generators=0 Angstrom
outfile2<< "Atomtypes="<<num_type<<" Charge="<<q_charge<<" Generators=0 Angstrom"<<endl;
//fprintf(outfile2,"Atomtypes=%i Charge=%.6f Generators=0 Angstrom",num_type,q_charge);
//print individual atom directives
int chk_tot=0;
for(int ll=1;ll<=num_type; ll++){
//get charge of atom nuc. 
int a_ch = get_charge(atom_type[ll]);
if(a_ch == 0){
outfile2.close();
return 1;
}
outfile2 << "Charge="<<a_ch<<" Atoms="<<atom_cnt[ll]<<" Basis="<<basis<<endl;
int chk_num=0;
for(int mm=1;mm<=NUM; mm++){ //find all atoms that match atom type
if(atom[mm].type == atom_type[ll]){ //matching atom found, print line
chk_tot++;
chk_num++;
outfile2<<atom[mm].type<<"    \t" \
<<fixed << setprecision(6) \
<<atom[mm].posx<<" \t"<<atom[mm].posy<<" \t"<<atom[mm].posz \
<<" \t 0 \t"<<chk_tot<<endl;
}
}
if(chk_num != atom_cnt[ll]){
cerr << "Error miss match in number of atoms of type "<<atom_type[ll]<<endl;
outfile2.close();
return 1;
}
}
if(chk_tot!=NUM){
cerr << "Error miss match in number of atoms found "<<NUM<<" vs output "<<chk_tot<<endl;
outfile2.close();
return 1;
}
cout <<NUM<<" atoms written to "<<OUTFILE<<endl;
outfile2.close();

return 0;	
}

int addatom_mol(int CNT1, string LINE1){
if(LINE1.length() < 40) return 1; //error in file or line
        string a2 = LINE1.substr(31,2);
        //a2.Trim();
        atom[CNT1].type = a2;
        atom[CNT1].posx = StoF(LINE1.substr(0,11));   
        atom[CNT1].posy = StoF(LINE1.substr(12,9));
        atom[CNT1].posz = StoF(LINE1.substr(22,9));
//cout << "atom added, type " << atom[CNT1].type << " pol "<< atom[CNT1].bondPOL << endl;
return 0;
}
int addatom_pdb(int CNT1, string LINE1){
if(LINE1.substr(0,3)=="END") return 1; 
if(LINE1.substr(0,4)!="ATOM") return 2;//????
if(LINE1.length() < 60) return 2; //error in file or line
        string a2 = LINE1.substr(13,1);
if(p2_flag==false){
string a3 = LINE1.substr(14,1);
if(a3=="1"){}      else if(a3=="2"){}
else if(a3=="3"){} else if(a3=="4"){}
else if(a3=="5"){} else if(a3=="6"){}
else if(a3=="7"){} else if(a3=="8"){}
else if(a3=="9"){} else if(a3=="0"){}
else{ a2 += a3;  } 
}
//a2.Trim(); cout << "Type detected "<<a2<<endl;
        atom[CNT1].type = a2;
        atom[CNT1].posx = StoF(LINE1.substr(31,7));   
        atom[CNT1].posy = StoF(LINE1.substr(39,7));
        atom[CNT1].posz = StoF(LINE1.substr(47,7));
string check2=LINE1.substr(30,30);
double comp[4];
std::istringstream iss(check2);
iss >> comp[1]; 
iss >> comp[2];
iss >> comp[3];
if((atom[CNT1].posx - comp[1]) > 0.01){
cout <<" Warning check x pos atom "<<CNT1<<" got "<<comp[1]<<" and "<<atom[CNT1].posx<<endl;
}
if((atom[CNT1].posy - comp[2]) > 0.01){
cout <<" Warning check y pos atom "<<CNT1<<" got "<<comp[2]<<" and "<<atom[CNT1].posy<<endl;
}

if((atom[CNT1].posz - comp[3]) > 0.01){
cout <<" Warning check z pos atom "<<CNT1<<" got "<<comp[3]<<" and "<<atom[CNT1].posz<<endl;
}

//cout << "atom added, type " << atom[CNT1].type << " pol "<< atom[CNT1].bondPOL << endl;
//cout <<","<<LINE1.substr(31,7)<<","<<LINE1.substr(39,7)<<","<<LINE1.substr(47,7)<<endl;
return 0;
}
int addatom_gro(int CNT1, string LINE1){
if(LINE1.length() < 40) return 1; //error in file or line
        string a2 = LINE1.substr(11,4);
        //a2.Trim();
        atom[CNT1].type = a2;
        atom[CNT1].posx = 10* StoF(LINE1.substr(22,6));   // X 10 for nm to Angstroms
        atom[CNT1].posy = 10* StoF(LINE1.substr(30,6));
        atom[CNT1].posz = 10* StoF(LINE1.substr(38,6));
//cout << "atom added, type " << atom[CNT1].type << " pol "<< atom[CNT1].bondPOL << endl;
return 0;
}
int addatom_xyz(int CNT1, string LINE1){
    string FRAG;
stringstream ss(LINE1);
vector<string> tokens;
int cc=1;
while( ss >> FRAG){
   tokens.push_back(FRAG);
if(cc==1) atom[CNT1].type =FRAG; //atomtype is 1st 
else if(cc==2) atom[CNT1].posx =StoF(FRAG); //x coord 
else if(cc==3) atom[CNT1].posy =StoF(FRAG); //y coord
else if(cc==4) atom[CNT1].posz =StoF(FRAG); //z coord
cc++;
}
if(cc<5) return 1; //error in file or line
return 0;
}

double StoF(string CONV){

  float val ;
  stringstream ss (stringstream::in | stringstream::out);
  ss << CONV;
  ss >> val ;
  //cout << "float" << val <<endl;
  return val;
}
bool is_number(char c_a){ //check to see if a character is a number 0-9
if(c_a=='1') return true;
if(c_a=='2') return true;
if(c_a=='3') return true;
if(c_a=='4') return true;
if(c_a=='5') return true;
if(c_a=='6') return true;
if(c_a=='7') return true;
if(c_a=='8') return true;
if(c_a=='9') return true;
if(c_a=='0') return true;
return false;
}

int get_charge(string mol_sym){
if(mol_sym== "H") return  1 ;
if(mol_sym== "He") return 2  ;
if(mol_sym== "Li") return 3  ;
if(mol_sym== "Be") return 4  ;
if(mol_sym== "B") return  5 ;
if(mol_sym== "C") return  6 ;
if(mol_sym== "N") return  7 ;
if(mol_sym== "O") return  8 ;
if(mol_sym== "F") return  9 ;
if(mol_sym== "Ne") return 10  ;
if(mol_sym== "Na") return 11  ;
if(mol_sym== "Mg") return 12  ;
if(mol_sym== "Al") return 13  ;
if(mol_sym== "Si") return 14  ;
if(mol_sym== "P") return  15 ;
if(mol_sym== "S") return  16 ;
if(mol_sym== "Cl") return 17  ;
if(mol_sym== "Ar") return 18  ;
if(mol_sym== "K") return  19 ;
if(mol_sym== "Ca") return 20  ;
if(mol_sym== "Sc") return 21  ;
if(mol_sym== "Ti") return 22  ;
if(mol_sym== "V") return  23 ;
if(mol_sym== "Cr") return 24  ;
if(mol_sym== "Mn") return 25  ;
if(mol_sym== "Fe") return 26  ;
if(mol_sym== "Co") return 27  ;
if(mol_sym== "Ni") return 28  ;
if(mol_sym== "Cu") return 29  ;
if(mol_sym== "Zn") return 30  ;
if(mol_sym== "Ga") return 31  ;
if(mol_sym== "Ge") return 32  ;
if(mol_sym== "As") return 33  ;
if(mol_sym== "Se") return 34  ;
if(mol_sym== "Br") return 35  ;
if(mol_sym== "Kr") return 36  ;
if(mol_sym== "Rb") return 37  ;
if(mol_sym== "Sr") return 38  ;
if(mol_sym== "Y") return  39 ;
if(mol_sym== "Zr") return 40  ;
if(mol_sym== "Nb") return 41  ;
if(mol_sym== "Mo") return 42  ;
if(mol_sym== "Tc") return 43  ;
if(mol_sym== "Ru") return 44  ;
if(mol_sym== "Rh") return 45  ;
if(mol_sym== "Pd") return 46  ;
if(mol_sym== "Ag") return 47  ;
if(mol_sym== "Cd") return 48  ;
if(mol_sym== "In") return 49  ;
if(mol_sym== "Sn") return 50  ;
if(mol_sym== "Sb") return 51  ;
if(mol_sym== "Te") return 52  ;
if(mol_sym== "I") return  53 ;
if(mol_sym== "Xe") return 54  ;
if(mol_sym== "Cs") return 55  ;
if(mol_sym== "Ba") return 56  ;
if(mol_sym== "La") return 57  ;
if(mol_sym== "Ce") return 58  ;
if(mol_sym== "Pr") return 59  ;
if(mol_sym== "Nd") return 60  ;
if(mol_sym== "Pm") return 61  ;
if(mol_sym== "Sm") return 62  ;
if(mol_sym== "Eu") return 63  ;
if(mol_sym== "Gd") return 64  ;
if(mol_sym== "Tb") return 65  ;
if(mol_sym== "Dy") return 66  ;
if(mol_sym== "Ho") return 67  ;
if(mol_sym== "Er") return 68  ;
if(mol_sym== "Tm") return 69  ;
if(mol_sym== "Yb") return 70  ;
if(mol_sym== "Lu") return 71  ;
if(mol_sym== "Hf") return 72  ;
if(mol_sym== "Ta") return 73  ;
if(mol_sym== "W") return  74 ;
if(mol_sym== "Re") return 75  ;
if(mol_sym== "Os") return 76  ;
if(mol_sym== "Ir") return 77  ;
if(mol_sym== "Pt") return 78  ;
if(mol_sym== "Au") return 79  ;
if(mol_sym== "Hg") return 80  ;
if(mol_sym== "Tl") return 81  ;
if(mol_sym== "Pb") return 82  ;
if(mol_sym== "Bi") return 83  ;
if(mol_sym== "Po") return 84  ;
if(mol_sym== "At") return 85  ;
if(mol_sym== "Rn") return 86  ;
if(mol_sym== "Fr") return 87  ;
if(mol_sym== "Ra") return 88  ;
if(mol_sym== "Ac") return 89  ;
if(mol_sym== "Th") return 90  ;
if(mol_sym== "Pa") return 91  ;
if(mol_sym== "U") return  92 ;	
cerr <<" undef atom type "<<mol_sym<<endl;
return 0;
}
