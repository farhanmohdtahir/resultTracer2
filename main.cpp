#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <cstring>
#include <sstream>
#include<locale>
using namespace std;

struct rng{
    string totalReadRng,noSubs, subs;
    int noMapCnt, softClipCnt, subsSoftClipCnt ;//totalReadCnt, noSubsCnt, subsCnt, misMapCnt, trueMapCnt,  trueMapCntSubs, misMapCntSubs, noMapCntSubs;
    int multiMapCnt, readMultiMapCnt, singleMapCnt, trueReadMultiMapCnt, misReadMultiMapCnt, trueMultiMapCnt, misMultiMapCnt, trueSingleMapCnt, misSingleMapCnt ;
    int subsMultiMapCnt, subsReadMultiMapCnt, subsSingleMapCnt, subsTrueReadMultiMapCnt, subsMisReadMultiMapCnt, subsTrueMultiMapCnt, subsMisMultiMapCnt;
    int subsTrueSingleMapCnt, subsMisSingleMapCnt, totalReadSam, noSubsSam, subsSam ;    
};

int main(){
    std::locale loc;
    ifstream infile1, infile2, infile3, infile4, infile5, infile6;
    ofstream outfile;
    int  totalRef=0,totalGff=0, i=0, j=0, k=0, l=0, m=0, n=0;
    string in1="expected_normal.txt", in2="output_hsa_normal.sam", in3="hsa.gff3", out="result.txt", line, lineTemp, lineIden, gffPos, readHeadStr;
    int lineIdenFirst[20000], lineIdenSec[20000], fastqTotRead;

    infile1.open(in1.c_str());
        while(getline(infile1, line)){
            if(line[0]=='@') ++totalRef;
        }
    infile1.close();
    
    infile4.open(in3.c_str());
        while(getline(infile4, line)){
            if(line[0]!='#') ++totalGff;
        }
        infile4.close();
    // cout<<totalGff;
    
    string header[totalRef], totalRead[totalRef], noSubsMap[totalRef], subsMap[totalRef], subsBaseNum[totalRef], id[totalRef], subsLineStr;
    std::string::size_type sz;
    // int totalReadCnt[totalRef], noSubsCnt[totalRef], subsCnt[totalRef], softClipCnt[totalRef], noMapCnt[totalRef], misMapCnt[totalRef];
    rng a[totalRef],b[totalRef],c[totalRef],d[totalRef],e[totalRef],f[totalRef];
    int subsLine[totalRef][30];
    string gffID[totalGff], gffIDStr;
    int gffPosStrt[totalGff], gffPosEnd[totalGff];

    infile2.open(in1.c_str());
        while(getline(infile2, line)){
            if(line[0]=='@'){
                line.erase(find(line.begin(), line.end(), '@'));
                header[i]=line;
                for (j=0; j<line.length(); j++){
                    if (line[j]==';')break;
                    id[i]+=line[j];
                }
            }
            else if(line[0]=='#'){
                j=0;
                for(k=0; k<line.length(); k++){
                    if (k>0){
                        if(line[k]!='+'){
                            subsLineStr+=line[k];
                            stringstream(subsLineStr) >> subsLine[i][j];
                        }
                        else{
                            ++j;
                            subsLineStr="";
                        }
                    }
                }
                // cout<<subsLine[i][j];
            }
            else if (line[0]=='\0'){}
            else{  
                k=0;
                for (j=0; j<line.length(); j++){
                    if (line[j]=='+')++k;
                    else {
                        switch(k){
                            case 0:
                                a[i].totalReadRng+=line[j];
                                break;
                            case 1:
                                b[i].totalReadRng+=line[j];
                                break;
                            case 2:
                                c[i].totalReadRng+=line[j];
                                break;
                            case 3:
                                d[i].totalReadRng+=line[j];
                                break;
                            case 4:
                                e[i].totalReadRng+=line[j];
                                break;
                            case 5:
                                f[i].totalReadRng+=line[j];
                                break;
                            case 6:
                                a[i].noSubs+=line[j];
                                break;
                            case 7:
                                b[i].noSubs+=line[j];
                                break;
                            case 8:
                                c[i].noSubs+=line[j];
                                break;
                            case 9:
                                d[i].noSubs+=line[j];
                                break;
                            case 10:
                                e[i].noSubs+=line[j];
                                break;
                            case 11:
                                f[i].noSubs+=line[j];
                                break;
                           case 12:
                                a[i].subs+=line[j];
                                break;
                            case 13:
                                b[i].subs+=line[j];
                                break;
                            case 14:
                                c[i].subs+=line[j];
                                break;
                            case 15:
                                d[i].subs+=line[j];
                                break;
                            case 16:
                                e[i].subs+=line[j];
                                break;
                            case 17:
                                f[i].subs+=line[j];
                                break;                            
                        }
                    }
                    // cout<<k<<" ";
                }
                ++i;
            }
            subsLineStr="";
        }

        // int totSamLine;

        // string simRead[totalRef][100];  
        lineTemp=header[0];
        
        infile2.close();
// //----------------------------------------------------------------------------------------------------------------------------
        infile5.open(in3.c_str());
        i=0;
        while(getline(infile5, line)){
            if (line[0]!='#'){
                k=0;
                for(j=0; j<line.length(); j++){
                    if (line[j]=='\t'){
                        if (k==3){
                            gffPosStrt[i]=stoi(gffPos);
                            gffPos="";
                        }
                        if(k==4){
                            gffPosEnd[i]=stoi(gffPos);
                            gffPos="";
                        }
                        ++k;
                        }
                    else {
                        if (k==3) gffPos+=line[j];
                        if (k==4) gffPos+=line[j];
                        if (k==8) gffID[i]+=line[j];
                    }
                }
                ++i;
            }
        }
        infile5.close();
k=0;
        for (i=0; i<totalGff; i++){
            k=0;
            gffIDStr=gffID[i];
            l=0;
            while(gffIDStr[l]!=';'){
                ++l;
                }
            gffID[i]="";
            for(m=l+1; m<gffIDStr.length(); m++){
                    gffID[i]+=gffIDStr[m];
            }
            gffIDStr=gffID[i];
            gffID[i]="";
            for(m=5; m<gffIDStr.length(); m++){
                if (gffIDStr[m]=='-')++k;
                if(k==3||gffIDStr[m]=='"')break;
                gffID[i]+=tolower(gffIDStr[m],loc);
            }

        }

//--------------------------------------------------------------------------------------------------------------------------------------------------------
int totalSamLine[totalRef];
int samMax [totalRef]={0} ;
string  tempHead[2];

i=0;
m=0;
        infile3.open(in2.c_str());
        while(getline(infile3, line)){
            if (line[0]=='h'&&line[1]=='s'){
                for (j=0; j<line.length(); j++){
                    if(m!=0){
                        if(line[j]=='\t')break;
                        tempHead[1]+=line[j];
                    }
                    else{
                        if(line[j]=='\t')break;
                        tempHead[0]+=line[j];
                    }
                }
                if(m!=0){         
                    if (tempHead[0]==tempHead[1]){
                        ++samMax[i];
                        tempHead[1]="";
                    }
                    else{
                        ++i;
                        tempHead[0]=tempHead[1];
                        tempHead[1]="";                    
                    }
                }
                m=1;
            }
        }
        infile3.close();
//----------------------------------------------------------------------------------------------------------------------
infile6.open(in2.c_str());
int max=0;
for (i=0; i<totalRef; i++){
    if(samMax[i]>max)max=samMax[i];
}
max=totalRef;                           //need to see later. Causing a segmentation fault

string ** readHead=new string *[totalRef];
int ** readLen=new int *[totalRef];
string ** softClip=new string *[totalRef];
string ** subsDet=new string *[totalRef];
string ** multiDet=new string *[totalRef];
int ** locat=new int *[totalRef];
int ** mulNum=new int *[totalRef];
string locatStr, readLenStr, multiDetStr, subsDetStr, softClipStr, mulNumStr, mulNumStrTemp;
cout<<totalRef<<" "<<max<<endl;
for (i = 0; i < max; i++) {
      readHead[i] = new string[max];
      readLen[i] = new int[max];
      softClip[i] = new string[max];
      subsDet[i] = new string[max];
      multiDet[i] = new string[max];
      locat[i] = new int[max];
      mulNum[i] = new int [max];
   }
// readHead[3000][max];
        i=0;l=0,n=0;
        while(getline(infile6, line)){
                            // cout<<line<<endl;
            if (line[0]=='h'&&line[1]=='s'){
                k=0;
                if(l>samMax[i]){
                    ++i;
                    l=0;
                }
                // cout<<i<<" "<<l<<endl;
                for (j=0; j<line.length(); j++){
                    if (line[j]=='\t'){
                        if (k==0){
                            readHeadStr=readHead[i][l];
                            readHead[i][l]="";
                            n=0;
                            for(m=0;m<readHeadStr.length();m++){
                                if(readHeadStr[m]=='-')++n;
                                if(n==3||readHeadStr[m]=='-'&&readHeadStr[m+1]=='5'&&readHeadStr[m+2]=='p'||readHeadStr[m]=='-'&&readHeadStr[m+1]=='3'&&readHeadStr[m+2]=='p')break;
                                readHead[i][l]+=tolower(readHeadStr[m],loc);
                            }
                            // cout<<readHead[i][l]<<endl;
                            // cout<<i<<" "<<l<<" "<<readHeadStr<<endl;
                            readHeadStr="";
                        }
                        if (k==3){
                            locat[i][l]=stoi(locatStr);
                            locatStr="";
                        }
                        if (k==5){
                            if(strchr(softClipStr.c_str(), 'S')) softClip[i][l]="1";
                            else softClip[i][l]="0";
                            softClipStr="";
                        }
                        if(k==9){
                            readLen[i][l]=readLenStr.length();
                            readLenStr="";
                        }
                        if (k==15){
                            if (strchr(subsDetStr.c_str(), 'A')||strchr(subsDetStr.c_str(), 'T')||strchr(subsDetStr.c_str(), 'G')||
                            strchr(subsDetStr.c_str(), 'C')||strchr(subsDetStr.c_str(), 'N'))subsDet[i][l]="1";
                            else subsDet[i][l]="0";
                            subsDetStr="";
                        }
                        if(k==16){
                            if (multiDetStr[0]=='C' && multiDetStr[1]=='C')multiDet[i][l]="2";
                            else if( multiDetStr[0]=='Z' && multiDetStr[1]=='S') multiDet[i][l]="1";
                            else multiDet[i][l]="0";
                            multiDetStr="";
                        }
                        if (k>16){
                        mulNumStrTemp="";
                            if (multiDet[i][l]=="1"){
                                if (k==17){
                                    for(m=5; m<mulNumStr.length();m++){
                                        mulNumStrTemp+=mulNumStr[m];
                                    }
                                    mulNum[i][l]=stoi(mulNumStrTemp);
                                }
                            }
                            else if (multiDet[i][l]=="2"){
                                if (k==19){
                                    for(m=5; m<mulNumStr.length();m++){
                                        mulNumStrTemp+=mulNumStr[m];
                                    }
                                    mulNum[i][l]=stoi(mulNumStrTemp);
                                }
                            }
                            mulNumStr="";
                        }
                        ++k;
                    }
                    else{
                        if (k==0)readHead[i][l]+=line[j];
                        if (k==3)locatStr+=line[j];
                        if (k==5)softClipStr+=line[j];
                        if (k==9)readLenStr+=line[j];
                        if (k>9){
                            if (k==15) subsDetStr+=line[j];
                            if (k>15){
                                if (k==16) multiDetStr+=line[j];
                                if (k>16){
                                    if (multiDet[i][l]=="1"){
                                        if (k==17)mulNumStr+=line[j];
                                    }
                                    else if (multiDet[i][l]=="2"){
                                        if (k==19)mulNumStr+=line[j];
                                    }
                                }
                            }
                        }
                    }
                }
                ++l;
            }
        }

        infile6.close();
//---------------------------------------------------------------------------------------------------------------------------
l=0;
k=0;
        for (i=0; i<totalRef; i++){
            for (j=0; j<=samMax[i]; j++){
                if (readLen[i][j]<=17){
                    ++a[l].totalReadSam;
                    if(subsDet[i][j]=="0"){
                        ++a[l].noSubsSam;
                        if(multiDet[i][j]=="1"||multiDet[i][j]=="2"){
                            ++a[l].readMultiMapCnt;
                            if(multiDet[i][j]=="1")++a[l].multiMapCnt;
                            for (k=0; k<totalGff; k++){
                                if (gffID[k]==readHead[i][j]){
                                    if(locat[i][j]>=gffPosStrt[k]&&locat[i][j]<=gffPosEnd[k]){
                                        if (multiDet[i][j]=="2"){
                                            ++a[l].trueReadMultiMapCnt;
                                            ++n;
                                        }
                                        if(multiDet[i][j]=="1"){
                                            ++a[l].trueReadMultiMapCnt;
                                            ++n;
                                            if(n==mulNum[i][j])++a[l].trueMultiMapCnt;
                                            n=0;
                                            }
                                        break;
                                    }
                                }
                            }
                        }
                        else if (multiDet[i][j]=="0"){
                            ++a[l].singleMapCnt;
                            for (k=0; k<totalGff; k++){
                                if (gffID[k]==readHead[i][j]){
                                    if(locat[i][j]>=gffPosStrt[k]&&locat[i][j]<=gffPosEnd[k]){
                                        ++a[l].trueSingleMapCnt;
                                        break;
                                    }
                                }
                            }
                        }
                        if(softClip[i][j]=="1")++a[l].softClipCnt;
                    }

                    else if(subsDet[i][j]=="1"){
                        ++a[l].subsSam;
                        if(multiDet[i][j]=="1"||multiDet[i][j]=="2"){
                            ++a[l].subsReadMultiMapCnt;
                            if(multiDet[i][j]=="1")++a[l].subsMultiMapCnt;
                            for (k=0; k<totalGff; k++){
                                if (gffID[k]==readHead[i][j]){
                                    if(locat[i][j]>=gffPosStrt[k]&&locat[i][j]<=gffPosEnd[k]){
                                        if (multiDet[i][j]=="2"){
                                            ++a[l].subsTrueReadMultiMapCnt;
                                            ++n;
                                        }
                                        if(multiDet[i][j]=="1"){
                                            ++a[l].subsTrueReadMultiMapCnt;
                                            ++n;
                                            if(n==mulNum[i][j])++a[l].subsTrueMultiMapCnt;
                                            n=0;
                                            }
                                        break;
                                    }
                                }
                            }
                        }
                        else if (multiDet[i][j]=="0"){
                            ++a[l].subsSingleMapCnt;
                            for (k=0; k<totalGff; k++){
                                if (gffID[k]==readHead[i][j]){
                                    if(locat[i][j]>=gffPosStrt[k]&&locat[i][j]<=gffPosEnd[k]){
                                        ++a[l].subsTrueSingleMapCnt;
                                        break;
                                    }
                                }
                            }
                        }
                        if(softClip[i][j]=="1")++a[l].subsSoftClipCnt;
                    }

                    else if(subsDet[i][j]==""){
                        ++a[l].noMapCnt;
                    }
                }

                if (readLen[i][j]==18){
                    ++b[l].totalReadSam;
                    if(subsDet[i][j]=="0"){
                        ++b[l].noSubsSam;
                        if(multiDet[i][j]=="1"||multiDet[i][j]=="2"){
                            ++b[l].readMultiMapCnt;
                            if(multiDet[i][j]=="1")++b[l].multiMapCnt;
                            for (k=0; k<totalGff; k++){
                                if (gffID[k]==readHead[i][j]){
                                    // cout<<i<<"\t"<<j<<"\t"<<gffPosStrt[k]<<"\t"<<locat[i][j]<<"\t"<<gffPosEnd[k]<<endl;
                                    if(locat[i][j]>=gffPosStrt[k]&&locat[i][j]<=gffPosEnd[k]){
                                        if (multiDet[i][j]=="2"){
                                            ++b[l].trueReadMultiMapCnt;
                                            ++n;
                                        }
                                        if(multiDet[i][j]=="1"){
                                            ++b[l].trueReadMultiMapCnt;
                                            ++n;
                                            if(n==mulNum[i][j])++b[l].trueMultiMapCnt;
                                            n=0;
                                            }
                                        break;
                                    }
                                }
                            }
                        }
                        else if (multiDet[i][j]=="0"){
                            ++b[l].singleMapCnt;
                            for (k=0; k<totalGff; k++){
                                if (gffID[k]==readHead[i][j]){
                                    if(locat[i][j]>=gffPosStrt[k]&&locat[i][j]<=gffPosEnd[k]){
                                        ++b[l].trueSingleMapCnt;
                                        break;
                                    }
                                }
                            }
                        }
                        if(softClip[i][j]=="1")++b[l].softClipCnt;
                    }

                    else if(subsDet[i][j]=="1"){
                        ++b[l].subsSam;
                        if(multiDet[i][j]=="1"||multiDet[i][j]=="2"){
                            ++b[l].subsReadMultiMapCnt;
                            if(multiDet[i][j]=="1")++b[l].subsMultiMapCnt;
                            for (k=0; k<totalGff; k++){
                                if (gffID[k]==readHead[i][j]){
                                    if(locat[i][j]>=gffPosStrt[k]&&locat[i][j]<=gffPosEnd[k]){
                                        if (multiDet[i][j]=="2"){
                                            ++b[l].subsTrueReadMultiMapCnt;
                                            ++n;
                                        }
                                        if(multiDet[i][j]=="1"){
                                            ++b[l].subsTrueReadMultiMapCnt;
                                            ++n;
                                            if(n==mulNum[i][j])++b[l].subsTrueMultiMapCnt;
                                            n=0;
                                            }
                                        break;
                                    }
                                }
                            }
                        }
                        else if (multiDet[i][j]=="0"){
                            ++b[l].subsSingleMapCnt;
                            for (k=0; k<totalGff; k++){
                                if (gffID[k]==readHead[i][j]){
                                    if(locat[i][j]>=gffPosStrt[k]&&locat[i][j]<=gffPosEnd[k]){
                                        ++b[l].subsTrueSingleMapCnt;
                                        break;
                                    }
                                }
                            }
                        }
                        if(softClip[i][j]=="1")++b[l].subsSoftClipCnt;
                    }

                    else if(subsDet[i][j]==""){
                        ++b[l].noMapCnt;
                    }
                }

                if (readLen[i][j]==19){
                    ++c[l].totalReadSam;
                    if(subsDet[i][j]=="0"){
                        ++c[l].noSubsSam;
                        if(multiDet[i][j]=="1"||multiDet[i][j]=="2"){
                            ++c[l].readMultiMapCnt;
                            if(multiDet[i][j]=="1")++c[l].multiMapCnt;
                            for (k=0; k<totalGff; k++){
                                if (gffID[k]==readHead[i][j]){
                                    if(locat[i][j]>=gffPosStrt[k]&&locat[i][j]<=gffPosEnd[k]){
                                        if (multiDet[i][j]=="2"){
                                            ++c[l].trueReadMultiMapCnt;
                                            ++n;
                                        }
                                        if(multiDet[i][j]=="1"){
                                            ++c[l].trueReadMultiMapCnt;
                                            ++n;
                                            if(n==mulNum[i][j])++c[l].trueMultiMapCnt;
                                            n=0;
                                            }
                                        break;
                                    }
                                }
                            }
                        }
                        else if (multiDet[i][j]=="0"){
                            ++c[l].singleMapCnt;
                            for (k=0; k<totalGff; k++){
                                if (gffID[k]==readHead[i][j]){
                                    if(locat[i][j]>=gffPosStrt[k]&&locat[i][j]<=gffPosEnd[k]){
                                        ++c[l].trueSingleMapCnt;
                                        break;
                                    }
                                }
                            }
                        }
                        if(softClip[i][j]=="1")++c[l].softClipCnt;
                    }

                    else if(subsDet[i][j]=="1"){
                        ++c[l].subsSam;
                        if(multiDet[i][j]=="1"||multiDet[i][j]=="2"){
                            ++c[l].subsReadMultiMapCnt;
                            if(multiDet[i][j]=="1")++c[l].subsMultiMapCnt;
                            for (k=0; k<totalGff; k++){
                                if (gffID[k]==readHead[i][j]){
                                    if(locat[i][j]>=gffPosStrt[k]&&locat[i][j]<=gffPosEnd[k]){
                                        if (multiDet[i][j]=="2"){
                                            ++c[l].subsTrueReadMultiMapCnt;
                                            ++n;
                                        }
                                        if(multiDet[i][j]=="1"){
                                            ++c[l].subsTrueReadMultiMapCnt;
                                            ++n;
                                            if(n==mulNum[i][j])++c[l].subsTrueMultiMapCnt;
                                            n=0;
                                            }
                                        break;
                                    }
                                }
                            }
                        }
                        else if (multiDet[i][j]=="0"){
                            ++c[l].subsSingleMapCnt;
                            for (k=0; k<totalGff; k++){
                                if (gffID[k]==readHead[i][j]){
                                    if(locat[i][j]>=gffPosStrt[k]&&locat[i][j]<=gffPosEnd[k]){
                                        ++c[l].subsTrueSingleMapCnt;
                                        break;
                                    }
                                }
                            }
                        }
                        if(softClip[i][j]=="1")++c[l].subsSoftClipCnt;
                    }

                    else if(subsDet[i][j]==""){
                        ++c[l].noMapCnt;
                    }
                }

                if (readLen[i][j]==20){
                    ++d[l].totalReadSam;
                    if(subsDet[i][j]=="0"){
                        ++d[l].noSubsSam;
                        if(multiDet[i][j]=="1"||multiDet[i][j]=="2"){
                            ++d[l].readMultiMapCnt;
                            if(multiDet[i][j]=="1")++d[l].multiMapCnt;
                            for (k=0; k<totalGff; k++){
                                if (gffID[k]==readHead[i][j]){
                                    if(locat[i][j]>=gffPosStrt[k]&&locat[i][j]<=gffPosEnd[k]){
                                        if (multiDet[i][j]=="2"){
                                            ++d[l].trueReadMultiMapCnt;
                                            ++n;
                                        }
                                        if(multiDet[i][j]=="1"){
                                            ++d[l].trueReadMultiMapCnt;
                                            ++n;
                                            if(n==mulNum[i][j])++d[l].trueMultiMapCnt;
                                            n=0;
                                            }
                                        break;
                                    }
                                }
                            }
                        }
                        else if (multiDet[i][j]=="0"){
                            ++d[l].singleMapCnt;
                            for (k=0; k<totalGff; k++){
                                if (gffID[k]==readHead[i][j]){
                                    if(locat[i][j]>=gffPosStrt[k]&&locat[i][j]<=gffPosEnd[k]){
                                        ++d[l].trueSingleMapCnt;
                                        break;
                                    }
                                }
                            }
                        }
                        if(softClip[i][j]=="1")++d[l].softClipCnt;
                    }

                    else if(subsDet[i][j]=="1"){
                        ++d[l].subsSam;
                        if(multiDet[i][j]=="1"||multiDet[i][j]=="2"){
                            ++d[l].subsReadMultiMapCnt;
                            if(multiDet[i][j]=="1")++d[l].subsMultiMapCnt;
                            for (k=0; k<totalGff; k++){
                                if (gffID[k]==readHead[i][j]){
                                    if(locat[i][j]>=gffPosStrt[k]&&locat[i][j]<=gffPosEnd[k]){
                                        if (multiDet[i][j]=="2"){
                                            ++d[l].subsTrueReadMultiMapCnt;
                                            ++n;
                                        }
                                        if(multiDet[i][j]=="1"){
                                            ++d[l].subsTrueReadMultiMapCnt;
                                            ++n;
                                            if(n==mulNum[i][j])++d[l].subsTrueMultiMapCnt;
                                            n=0;
                                            }
                                        break;
                                    }
                                }
                            }
                        }
                        else if (multiDet[i][j]=="0"){
                            ++d[l].subsSingleMapCnt;
                            for (k=0; k<totalGff; k++){
                                if (gffID[k]==readHead[i][j]){
                                    if(locat[i][j]>=gffPosStrt[k]&&locat[i][j]<=gffPosEnd[k]){
                                        ++d[l].subsTrueSingleMapCnt;
                                        break;
                                    }
                                }
                            }
                        }
                        if(softClip[i][j]=="1")++d[l].subsSoftClipCnt;
                    }

                    else if(subsDet[i][j]==""){
                        ++d[l].noMapCnt;
                    }
                }

                if (readLen[i][j]==21){
                    ++e[l].totalReadSam;
                    if(subsDet[i][j]=="0"){
                        ++e[l].noSubsSam;
                        if(multiDet[i][j]=="1"||multiDet[i][j]=="2"){
                            ++e[l].readMultiMapCnt;
                            if(multiDet[i][j]=="1")++e[l].multiMapCnt;
                            for (k=0; k<totalGff; k++){
                                if (gffID[k]==readHead[i][j]){
                                    if(locat[i][j]>=gffPosStrt[k]&&locat[i][j]<=gffPosEnd[k]){
                                        if (multiDet[i][j]=="2"){
                                            ++e[l].trueReadMultiMapCnt;
                                            ++n;
                                        }
                                        if(multiDet[i][j]=="1"){
                                            ++e[l].trueReadMultiMapCnt;
                                            ++n;
                                            if(n==mulNum[i][j])++e[l].trueMultiMapCnt;
                                            n=0;
                                            }
                                        break;
                                    }
                                }
                            }
                        }
                        else if (multiDet[i][j]=="0"){
                            ++e[l].singleMapCnt;
                            for (k=0; k<totalGff; k++){
                                if (gffID[k]==readHead[i][j]){
                                    if(locat[i][j]>=gffPosStrt[k]&&locat[i][j]<=gffPosEnd[k]){
                                        ++e[l].trueSingleMapCnt;
                                        break;
                                    }
                                }
                            }
                        }
                        if(softClip[i][j]=="1")++e[l].softClipCnt;
                    }

                    else if(subsDet[i][j]=="1"){
                        ++e[l].subsSam;
                        if(multiDet[i][j]=="1"||multiDet[i][j]=="2"){
                            ++e[l].subsReadMultiMapCnt;
                            if(multiDet[i][j]=="1")++e[l].subsMultiMapCnt;
                            for (k=0; k<totalGff; k++){
                                if (gffID[k]==readHead[i][j]){
                                    if(locat[i][j]>=gffPosStrt[k]&&locat[i][j]<=gffPosEnd[k]){
                                        if (multiDet[i][j]=="2"){
                                            ++e[l].subsTrueReadMultiMapCnt;
                                            ++n;
                                        }
                                        if(multiDet[i][j]=="1"){
                                            ++e[l].subsTrueReadMultiMapCnt;
                                            ++n;
                                            if(n==mulNum[i][j])++e[l].subsTrueMultiMapCnt;
                                            n=0;
                                            }
                                        break;
                                    }
                                }
                            }
                        }
                        else if (multiDet[i][j]=="0"){
                            ++e[l].subsSingleMapCnt;
                            for (k=0; k<totalGff; k++){
                                if (gffID[k]==readHead[i][j]){
                                    if(locat[i][j]>=gffPosStrt[k]&&locat[i][j]<=gffPosEnd[k]){
                                        ++e[l].subsTrueSingleMapCnt;
                                        break;
                                    }
                                }
                            }
                        }
                        if(softClip[i][j]=="1")++e[l].subsSoftClipCnt;
                    }

                    else if(subsDet[i][j]==""){
                        ++e[l].noMapCnt;
                    }
                }

                if (readLen[i][j]==22){
                    ++f[l].totalReadSam;
                    if(subsDet[i][j]=="0"){
                        ++f[l].noSubsSam;
                        if(multiDet[i][j]=="1"||multiDet[i][j]=="2"){
                            ++f[l].readMultiMapCnt;
                            if(multiDet[i][j]=="1")++f[l].multiMapCnt;
                            for (k=0; k<totalGff; k++){
                                if (gffID[k]==readHead[i][j]){
                                    if(locat[i][j]>=gffPosStrt[k]&&locat[i][j]<=gffPosEnd[k]){
                                        if (multiDet[i][j]=="2"){
                                            ++f[l].trueReadMultiMapCnt;
                                            ++n;
                                        }
                                        if(multiDet[i][j]=="1"){
                                            ++f[l].trueReadMultiMapCnt;
                                            ++n;
                                            if(n==mulNum[i][j])++f[l].trueMultiMapCnt;
                                            n=0;
                                            }
                                        break;
                                    }
                                }
                            }
                        }
                        else if (multiDet[i][j]=="0"){
                            ++f[l].singleMapCnt;
                            for (k=0; k<totalGff; k++){
                                if (gffID[k]==readHead[i][j]){
                                    if(locat[i][j]>=gffPosStrt[k]&&locat[i][j]<=gffPosEnd[k]){
                                        ++f[l].trueSingleMapCnt;
                                        break;
                                    }
                                }
                            }
                        }
                        if(softClip[i][j]=="1")++f[l].softClipCnt;
                    }

                    else if(subsDet[i][j]=="1"){
                        ++f[l].subsSam;
                        if(multiDet[i][j]=="1"||multiDet[i][j]=="2"){
                            ++f[l].subsReadMultiMapCnt;
                            if(multiDet[i][j]=="1")++f[l].subsMultiMapCnt;
                            for (k=0; k<totalGff; k++){
                                if (gffID[k]==readHead[i][j]){
                                    if(locat[i][j]>=gffPosStrt[k]&&locat[i][j]<=gffPosEnd[k]){
                                        if (multiDet[i][j]=="2"){
                                            ++f[l].subsTrueReadMultiMapCnt;
                                            ++n;
                                        }
                                        if(multiDet[i][j]=="1"){
                                            ++f[l].subsTrueReadMultiMapCnt;
                                            ++n;
                                            if(n==mulNum[i][j])++f[l].subsTrueMultiMapCnt;
                                            n=0;
                                            }
                                        break;
                                    }
                                }
                            }
                        }
                        else if (multiDet[i][j]=="0"){
                            ++f[l].subsSingleMapCnt;
                            for (k=0; k<totalGff; k++){
                                if (gffID[k]==readHead[i][j]){
                                    if(locat[i][j]>=gffPosStrt[k]&&locat[i][j]<=gffPosEnd[k]){
                                        ++f[l].subsTrueSingleMapCnt;
                                        break;
                                    }
                                }
                            }
                        }
                        if(softClip[i][j]=="1")++f[l].subsSoftClipCnt;
                    }

                    else if(subsDet[i][j]==""){
                        ++f[l].noMapCnt;
                    }
                }
            }
            ++l;
        }


    for (i=0;i<totalRef; i++){
        a[i].misMultiMapCnt=a[i].multiMapCnt-a[i].trueMultiMapCnt;
        b[i].misMultiMapCnt=b[i].multiMapCnt-b[i].trueMultiMapCnt;
        c[i].misMultiMapCnt=c[i].multiMapCnt-c[i].trueMultiMapCnt;
        d[i].misMultiMapCnt=d[i].multiMapCnt-d[i].trueMultiMapCnt;
        e[i].misMultiMapCnt=e[i].multiMapCnt-e[i].trueMultiMapCnt;
        f[i].misMultiMapCnt=f[i].multiMapCnt-f[i].trueMultiMapCnt;

        a[i].misReadMultiMapCnt=a[i].readMultiMapCnt-a[i].trueReadMultiMapCnt;
        b[i].misReadMultiMapCnt=b[i].readMultiMapCnt-b[i].trueReadMultiMapCnt;
        c[i].misReadMultiMapCnt=c[i].readMultiMapCnt-c[i].trueReadMultiMapCnt;
        d[i].misReadMultiMapCnt=d[i].readMultiMapCnt-d[i].trueReadMultiMapCnt;
        e[i].misReadMultiMapCnt=e[i].readMultiMapCnt-e[i].trueReadMultiMapCnt;
        f[i].misReadMultiMapCnt=f[i].readMultiMapCnt-f[i].trueReadMultiMapCnt;

        a[i].misSingleMapCnt=a[i].singleMapCnt-a[i].trueSingleMapCnt;
        b[i].misSingleMapCnt=b[i].singleMapCnt-b[i].trueSingleMapCnt;
        c[i].misSingleMapCnt=c[i].singleMapCnt-c[i].trueSingleMapCnt;
        d[i].misSingleMapCnt=d[i].singleMapCnt-d[i].trueSingleMapCnt;
        e[i].misSingleMapCnt=e[i].singleMapCnt-e[i].trueSingleMapCnt;
        f[i].misSingleMapCnt=f[i].singleMapCnt-f[i].trueSingleMapCnt;
    }

    
    for (i=0; i<totalRef; i++){
        cout<<header[i]<<endl<<"--------------------"<<endl;
        cout<<"\t\t17\t\t\t18\t\t\t19\t\t\t20\t\t\t21\t\t\t22"<<endl;
        cout<<"\t\tFASTQseq\tsamSeq\tFASTQseq\tsamSeq\tFASTQseq\tsamSeq\tFASTQseq\tsamSeq\tFASTQseq\tsamSeq\tFASTQseq\tsamSeq"<<endl;
        cout<<"total:\t\t"<<a[i].totalReadRng<<"\t\t"<<a[i].totalReadSam<<"\t"<<b[i].totalReadRng<<"\t\t"<<b[i].totalReadSam<<"\t"<<c[i].totalReadRng<<"\t\t"<<c[i].totalReadSam<<"\t";
        cout<<d[i].totalReadRng<<"\t\t"<<d[i].totalReadSam<<"\t"<<e[i].totalReadRng<<"\t\t"<<e[i].totalReadSam<<"\t"<<f[i].totalReadRng<<"\t\t"<<f[i].totalReadSam<<endl;

        cout<<"NoSubs:\t\t"<<a[i].noSubs<<"\t\t"<<a[i].noSubsSam<<"\t"<<b[i].noSubs<<"\t\t"<<b[i].noSubsSam<<"\t"<<c[i].noSubs<<"\t\t"<<c[i].noSubsSam<<"\t";
        cout<<d[i].noSubs<<"\t\t"<<d[i].noSubsSam<<"\t"<<e[i].noSubs<<"\t\t"<<e[i].noSubsSam<<"\t"<<f[i].noSubs<<"\t\t"<<f[i].noSubsSam<<endl;       

        cout<<"MULTIMAP:\t"<<a[i].multiMapCnt<<"\t\t"<<a[i].readMultiMapCnt<<"\t"<<b[i].multiMapCnt<<"\t\t"<<b[i].readMultiMapCnt<<"\t"<<c[i].multiMapCnt<<"\t\t"<<c[i].readMultiMapCnt<<"\t";
        cout<<d[i].multiMapCnt<<"\t\t"<<d[i].readMultiMapCnt<<"\t"<<e[i].multiMapCnt<<"\t\t"<<e[i].readMultiMapCnt<<"\t"<<f[i].multiMapCnt<<"\t\t"<<f[i].readMultiMapCnt<<endl;

        cout<<"correct:\t"<<a[i].trueMultiMapCnt<<"\t\t"<<a[i].trueReadMultiMapCnt<<"\t"<<b[i].trueMultiMapCnt<<"\t\t"<<b[i].trueReadMultiMapCnt<<"\t"<<c[i].trueMultiMapCnt<<"\t\t"<<c[i].trueReadMultiMapCnt<<"\t";
        cout<<d[i].trueMultiMapCnt<<"\t\t"<<d[i].trueReadMultiMapCnt<<"\t"<<e[i].trueMultiMapCnt<<"\t\t"<<e[i].trueReadMultiMapCnt<<"\t"<<f[i].trueMultiMapCnt<<"\t\t"<<f[i].trueReadMultiMapCnt<<endl;

        cout<<"incorrect:\t"<<a[i].misMultiMapCnt<<"\t\t"<<a[i].misReadMultiMapCnt<<"\t"<<b[i].misMultiMapCnt<<"\t\t"<<b[i].misReadMultiMapCnt<<"\t"<<c[i].misMultiMapCnt<<"\t\t"<<c[i].misReadMultiMapCnt<<"\t";
        cout<<d[i].misMultiMapCnt<<"\t\t"<<d[i].misReadMultiMapCnt<<"\t"<<e[i].misMultiMapCnt<<"\t\t"<<e[i].misReadMultiMapCnt<<"\t"<<f[i].misMultiMapCnt<<"\t\t"<<f[i].misReadMultiMapCnt<<endl;                          

        cout<<"SINGLEMAP:\t"<<a[i].singleMapCnt<<"\t\t"<<a[i].singleMapCnt<<"\t"<<b[i].singleMapCnt<<"\t\t"<<b[i].singleMapCnt<<"\t"<<c[i].singleMapCnt<<"\t\t"<<c[i].singleMapCnt<<"\t";
        cout<<d[i].singleMapCnt<<"\t\t"<<d[i].singleMapCnt<<"\t"<<e[i].singleMapCnt<<"\t\t"<<e[i].singleMapCnt<<"\t"<<f[i].singleMapCnt<<"\t\t"<<f[i].singleMapCnt<<endl;

        cout<<"correct:\t"<<a[i].trueSingleMapCnt<<"\t\t"<<a[i].trueSingleMapCnt<<"\t"<<b[i].trueSingleMapCnt<<"\t\t"<<b[i].trueSingleMapCnt<<"\t"<<c[i].trueSingleMapCnt<<"\t\t"<<c[i].trueSingleMapCnt<<"\t";
        cout<<d[i].trueSingleMapCnt<<"\t\t"<<d[i].trueSingleMapCnt<<"\t"<<e[i].trueSingleMapCnt<<"\t\t"<<e[i].trueSingleMapCnt<<"\t"<<f[i].trueSingleMapCnt<<"\t\t"<<f[i].trueSingleMapCnt<<endl;

        cout<<"incorrect:\t"<<a[i].misSingleMapCnt<<"\t\t"<<a[i].misSingleMapCnt<<"\t"<<b[i].misSingleMapCnt<<"\t\t"<<b[i].misSingleMapCnt<<"\t"<<c[i].misSingleMapCnt<<"\t\t"<<c[i].misSingleMapCnt<<"\t";
        cout<<d[i].misSingleMapCnt<<"\t\t"<<d[i].misSingleMapCnt<<"\t"<<e[i].misSingleMapCnt<<"\t\t"<<e[i].misSingleMapCnt<<"\t"<<f[i].misSingleMapCnt<<"\t\t"<<f[i].misSingleMapCnt<<endl;
//----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
        cout<<"SUBS:\t\t"<<a[i].subs<<"\t\t"<<a[i].subsSam<<"\t"<<b[i].subs<<"\t\t"<<b[i].subsSam<<"\t"<<c[i].subs<<"\t\t"<<c[i].subsSam<<"\t";
        cout<<d[i].subs<<"\t\t"<<d[i].subsSam<<"\t"<<e[i].subs<<"\t\t"<<e[i].subsSam<<"\t"<<f[i].subs<<"\t\t"<<f[i].subsSam<<endl;       

        cout<<"MULTIMAP:\t"<<a[i].subsMultiMapCnt<<"\t\t"<<a[i].subsReadMultiMapCnt<<"\t"<<b[i].subsMultiMapCnt<<"\t\t"<<b[i].subsReadMultiMapCnt<<"\t"<<c[i].subsMultiMapCnt<<"\t\t"<<c[i].subsReadMultiMapCnt<<"\t";
        cout<<d[i].subsMultiMapCnt<<"\t\t"<<d[i].subsReadMultiMapCnt<<"\t"<<e[i].subsMultiMapCnt<<"\t\t"<<e[i].subsReadMultiMapCnt<<"\t"<<f[i].subsMultiMapCnt<<"\t\t"<<f[i].subsReadMultiMapCnt<<endl;

        cout<<"correct:\t"<<a[i].subsTrueMultiMapCnt<<"\t\t"<<a[i].subsTrueReadMultiMapCnt<<"\t"<<b[i].subsTrueMultiMapCnt<<"\t\t"<<b[i].subsTrueReadMultiMapCnt<<"\t"<<c[i].subsTrueMultiMapCnt<<"\t\t"<<c[i].subsTrueReadMultiMapCnt<<"\t";
        cout<<d[i].subsTrueMultiMapCnt<<"\t\t"<<d[i].subsTrueReadMultiMapCnt<<"\t"<<e[i].subsTrueMultiMapCnt<<"\t\t"<<e[i].subsTrueReadMultiMapCnt<<"\t"<<f[i].subsTrueMultiMapCnt<<"\t\t"<<f[i].subsTrueReadMultiMapCnt<<endl;

        cout<<"incorrect:\t"<<a[i].subsMisMultiMapCnt<<"\t\t"<<a[i].subsMisReadMultiMapCnt<<"\t"<<b[i].subsMisMultiMapCnt<<"\t\t"<<b[i].subsMisReadMultiMapCnt<<"\t"<<c[i].subsMisMultiMapCnt<<"\t\t"<<c[i].subsMisReadMultiMapCnt<<"\t";
        cout<<d[i].subsMisMultiMapCnt<<"\t\t"<<d[i].subsMisReadMultiMapCnt<<"\t"<<e[i].subsMisMultiMapCnt<<"\t\t"<<e[i].subsMisReadMultiMapCnt<<"\t"<<f[i].subsMisMultiMapCnt<<"\t\t"<<f[i].subsMisReadMultiMapCnt<<endl;                          

        cout<<"SINGLEMAP:\t"<<a[i].subsSingleMapCnt<<"\t\t"<<a[i].subsSingleMapCnt<<"\t"<<b[i].subsSingleMapCnt<<"\t\t"<<b[i].subsSingleMapCnt<<"\t"<<c[i].subsSingleMapCnt<<"\t\t"<<c[i].subsSingleMapCnt<<"\t";
        cout<<d[i].subsSingleMapCnt<<"\t\t"<<d[i].subsSingleMapCnt<<"\t"<<e[i].subsSingleMapCnt<<"\t\t"<<e[i].subsSingleMapCnt<<"\t"<<f[i].subsSingleMapCnt<<"\t\t"<<f[i].subsSingleMapCnt<<endl;

        cout<<"correct:\t"<<a[i].subsTrueSingleMapCnt<<"\t\t"<<a[i].subsTrueSingleMapCnt<<"\t"<<b[i].subsTrueSingleMapCnt<<"\t\t"<<b[i].subsTrueSingleMapCnt<<"\t"<<c[i].subsTrueSingleMapCnt<<"\t\t"<<c[i].subsTrueSingleMapCnt<<"\t";
        cout<<d[i].subsTrueSingleMapCnt<<"\t\t"<<d[i].subsTrueSingleMapCnt<<"\t"<<e[i].subsTrueSingleMapCnt<<"\t\t"<<e[i].subsTrueSingleMapCnt<<"\t"<<f[i].subsTrueSingleMapCnt<<"\t\t"<<f[i].subsTrueSingleMapCnt<<endl;

        cout<<"incorrect:\t"<<a[i].subsMisSingleMapCnt<<"\t\t"<<a[i].subsMisSingleMapCnt<<"\t"<<b[i].subsMisSingleMapCnt<<"\t\t"<<b[i].subsMisSingleMapCnt<<"\t"<<c[i].subsMisSingleMapCnt<<"\t\t"<<c[i].subsMisSingleMapCnt<<"\t";
        cout<<d[i].subsMisSingleMapCnt<<"\t\t"<<d[i].subsMisSingleMapCnt<<"\t"<<e[i].subsMisSingleMapCnt<<"\t\t"<<e[i].subsMisSingleMapCnt<<"\t"<<f[i].subsMisSingleMapCnt<<"\t\t"<<f[i].subsMisSingleMapCnt<<endl;

        cout<<"NOTMAP:\t\t"<<a[i].noMapCnt<<"\t\t"<<a[i].noMapCnt<<"\t"<<b[i].noMapCnt<<"\t\t"<<b[i].noMapCnt<<"\t"<<c[i].noMapCnt<<"\t\t"<<c[i].noMapCnt<<"\t";
        cout<<d[i].noMapCnt<<"\t\t"<<d[i].noMapCnt<<"\t"<<e[i].noMapCnt<<"\t\t"<<e[i].noMapCnt<<"\t"<<f[i].noMapCnt<<"\t\t"<<f[i].noMapCnt<<endl<<endl;
        // cout<<"total:\t"<<a[i].multiMapCnt<<"\t\t"<<a[i].readMultiMapCnt<<"\t"<<b[i].multiMapCnt<<"\t\t"<<b[i].readMultiMapCnt<<"\t"<<c[i].multiMapCnt<<"\t\t"<<c[i].readMultiMapCnt<<"\t";
        // cout<<d[i].multiMapCnt<<"\t\t"<<d[i].readMultiMapCnt<<"\t"<<e[i].multiMapCnt<<"\t\t"<<e[i].readMultiMapCnt<<"\t"<<f[i].multiMapCnt<<"\t\t"<<f[i].readMultiMapCnt<<endl<<endl;
        // cout<<"Total Read:\t"<<totalRead[i]<<"\t|Total Read:\t\t"<<totalReadCnt[i]<<"\t|"<<endl;
        // cout<<"No Subs Read:\t"<<noSubsMap[i]<<"\t|No Subs Read Mapped:\t"<<noSubsCnt[i]<<"\t|"<<endl;
        // cout<<"Subs Read:\t"<<subsMap[i]<<"\t|SoftClip:\t\t"<<softClipCnt[i]<<"\t|"<<endl;
        // cout<<"\t\t\t|Subs Read Mapped:\t"<<softClipCnt[i]+subsCnt[i]<<"\t|"<<endl;
        // cout<<"Num Subs:\t"<<subsBaseNum[i]<<"\t|read not mapped:\t"<<noMapCnt[i]<<"\t|"<<endl;
        // cout<<"\t\t\t|incorrect mapped:\t"<<misMapCnt[i]<<"\t|"<<endl<<endl;
    }

// for (i=0; i<totalRef; i++){
//     for(j=0; j<samMax[i]; j++){
//         cout<<i<<"\t"<<j<<"\t"<<readHead[i][j]<<"\t"<<locat[i][j]<<"\t"<<softClip[i][j]<<"\t"<<readLen[i][j]<<"\t"<<subsDet[i][j]<<"\t"<<multiDet[i][j]<<"\t"<<mulNum[i][j]<<endl;
//     }
//     if(readHead[i][0]=="")break;    
//     cout<<endl;
// }

// for (i=0; i<totalGff; i++){
//     cout<<gffID[i]<<endl;
// }

return 0;
}
// for (i=0; i<totalRef; i++){
//     cout<<totalSamLine[i]<<endl;
// }



      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
        // for(i=0; i<totalRef; i++){
        //     for(j=0; j<totalSamLine[i]; j++){
        //         cout<<readHead[i][l]<<endl;
        //     }
        // }
 
    // for (i=0; i<totalRef; i++){cout<<i<<"\t"<<softClip[i][3]<<endl;}
//-----------------------------------------------------------------------------------------------------------------------------
        // j=0;
        // n=0;

        // for (i=0; i<totalRef; i++){
        //     totSamLine=stoi(a[i].totalReadRng)+stoi(b[i].totalReadRng)+stoi(c[i].totalReadRng)+stoi(d[i].totalReadRng)+stoi(e[i].totalReadRng)+stoi(f[i].totalReadRng);
        //     m=0;
        //     for (int l=0; l<totSamLine; l++){
                // if (readLen[i][l].length()==17){
                //                             // cout<<i<<" "<<l<<endl;
                //      if(readHead[i][l]==id[j]){
                //         ++a[j].totalReadCnt;
                //         if(l==subsLine[i][m]){
                //             if (strchr(subsDet[i][l].c_str(), 'A')||strchr(subsDet[i][l].c_str(), 'T')||strchr(subsDet[i][l].c_str(), 'G')||
                //             strchr(subsDet[i][l].c_str(), 'C')||strchr(subsDet[i][l].c_str(), 'N')) ++a[j].trueMapCntSubs;
                //             else if (i==lineIdenFirst[n]&&l==lineIdenSec[n]){
                //                 ++a[j].misMapCntSubs;
                //                 ++n;
                //             }                     
                //             else if (strchr(softClip[i][l].c_str(), 'S')) ++a[j].softClipCntSubs;
                //             else if (softClip[i][l]=="*") ++a[j].noMapCntSubs;
                //            // else ++a[j].trueMapCntSubs;  
                //             ++m;
                //         }
                //         else{
                //             if (strchr(subsDet[i][l].c_str(), 'A')||strchr(subsDet[i][l].c_str(), 'T')||strchr(subsDet[i][l].c_str(), 'G')||
                //             strchr(subsDet[i][l].c_str(), 'C')||strchr(subsDet[i][l].c_str(), 'N')) {}//++a[j].trueMapCnt;
                //             else if (i==lineIdenFirst[n]&&l==lineIdenSec[n]){
                //                 ++a[j].misMapCnt;
                //                 ++n;
                //             }                     
                //             else if (strchr(softClip[i][l].c_str(), 'S')) ++a[j].softClipCnt;//if(i==1){cout<<l<<","<<subsLine[i][m]<<" "<<i<<m<<"\t";}}
                //             else if (softClip[i][l]=="*") ++a[j].noMapCnt;
                //             else ++a[j].trueMapCnt;       
                //         }
                                    
                //     }
                // }

                // if (readLen[i][l].length()==18){
                //                             // cout<<i<<" "<<l<<endl;
                //      if(readHead[i][l]==id[j]){
                //         ++b[j].totalReadCnt;
                //         if(l==subsLine[i][m]){
                //             if (strchr(subsDet[i][l].c_str(), 'A')||strchr(subsDet[i][l].c_str(), 'T')||strchr(subsDet[i][l].c_str(), 'G')||
                //             strchr(subsDet[i][l].c_str(), 'C')||strchr(subsDet[i][l].c_str(), 'N')) ++b[j].trueMapCntSubs;
                //             else if (i==lineIdenFirst[n]&&l==lineIdenSec[n]){
                //                 ++b[j].misMapCntSubs;
                //                 ++n;
                //             }                     
                //             else if (strchr(softClip[i][l].c_str(), 'S')) ++b[j].softClipCntSubs;
                //             else if (softClip[i][l]=="*") ++b[j].noMapCntSubs;
                //             //else ++b[j].trueMapCntSubs;  
                //             ++m;
                //         }
                //         else{
                //             if (strchr(subsDet[i][l].c_str(), 'A')||strchr(subsDet[i][l].c_str(), 'T')||strchr(subsDet[i][l].c_str(), 'G')||
                //             strchr(subsDet[i][l].c_str(), 'C')||strchr(subsDet[i][l].c_str(), 'N')){} //++b[j].trueMapCnt;
                //             else if (i==lineIdenFirst[n]&&l==lineIdenSec[n]){
                //                 ++b[j].misMapCnt;
                //                 ++n;
                //             }                     
                //             else if (strchr(softClip[i][l].c_str(), 'S')) ++b[j].softClipCnt;
                //             else if (softClip[i][l]=="*") ++b[j].noMapCnt;
                //             else ++b[j].trueMapCnt;       
                //         }
                //     }
                // }

                // if (readLen[i][l].length()==19){
                //                             // cout<<i<<" "<<l<<endl;
                //      if(readHead[i][l]==id[j]){
                //         ++c[j].totalReadCnt;
                //         if(l==subsLine[i][m]){
                //             if (strchr(subsDet[i][l].c_str(), 'A')||strchr(subsDet[i][l].c_str(), 'T')||strchr(subsDet[i][l].c_str(), 'G')||
                //             strchr(subsDet[i][l].c_str(), 'C')||strchr(subsDet[i][l].c_str(), 'N')) ++c[j].trueMapCntSubs;
                //             else if (i==lineIdenFirst[n]&&l==lineIdenSec[n]){
                //                 ++c[j].misMapCntSubs;
                //                 ++n;
                //             }                     
                //             else if (strchr(softClip[i][l].c_str(), 'S')) ++c[j].softClipCntSubs;
                //             else if (softClip[i][l]=="*") ++c[j].noMapCntSubs;
                //             //else ++c[j].trueMapCntSubs;  
                //             ++m;
                //         }
                //         else{
                //             if (strchr(subsDet[i][l].c_str(), 'A')||strchr(subsDet[i][l].c_str(), 'T')||strchr(subsDet[i][l].c_str(), 'G')||
                //             strchr(subsDet[i][l].c_str(), 'C')||strchr(subsDet[i][l].c_str(), 'N')){}// ++c[j].trueMapCnt;
                //             else if (i==lineIdenFirst[n]&&l==lineIdenSec[n]){
                //                 ++c[j].misMapCnt;
                //                 ++n;
                //             }                     
                //             else if (strchr(softClip[i][l].c_str(), 'S')) ++c[j].softClipCnt;
                //             else if (softClip[i][l]=="*") ++c[j].noMapCnt;
                //             else ++c[j].trueMapCnt;       
                //         }
                //     }
                // }

                // if (readLen[i][l].length()==20){
                //                             // cout<<i<<" "<<l<<endl;
                //      if(readHead[i][l]==id[j]){
                //         ++d[j].totalReadCnt;
                //         if(l==subsLine[i][m]){
                //             if (strchr(subsDet[i][l].c_str(), 'A')||strchr(subsDet[i][l].c_str(), 'T')||strchr(subsDet[i][l].c_str(), 'G')||
                //             strchr(subsDet[i][l].c_str(), 'C')||strchr(subsDet[i][l].c_str(), 'N')) ++d[j].trueMapCntSubs;
                //             else if (i==lineIdenFirst[n]&&l==lineIdenSec[n]){
                //                 ++d[j].misMapCntSubs;
                //                 ++n;
                //             }                     
                //             else if (strchr(softClip[i][l].c_str(), 'S')) ++d[j].softClipCntSubs;
                //             else if (softClip[i][l]=="*") ++d[j].noMapCntSubs;
                //            // else ++d[j].trueMapCntSubs;  
                //             ++m;
                //         }
                //         else{
                //             if (strchr(subsDet[i][l].c_str(), 'A')||strchr(subsDet[i][l].c_str(), 'T')||strchr(subsDet[i][l].c_str(), 'G')||
                //             strchr(subsDet[i][l].c_str(), 'C')||strchr(subsDet[i][l].c_str(), 'N')) {}//++d[j].trueMapCnt;
                //             else if (i==lineIdenFirst[n]&&l==lineIdenSec[n]){
                //                 ++d[j].misMapCnt;
                //                 ++n;
                //             }                     
                //             else if (strchr(softClip[i][l].c_str(), 'S')) ++d[j].softClipCnt;
                //             else if (softClip[i][l]=="*") ++d[j].noMapCnt;
                //             else ++d[j].trueMapCnt;       
                //         }
                //     }
                // }

                // if (readLen[i][l].length()==21){
                //                             // cout<<i<<" "<<l<<endl;
                //      if(readHead[i][l]==id[j]){
                //         ++e[j].totalReadCnt;
                //         if(l==subsLine[i][m]){
                //             if (strchr(subsDet[i][l].c_str(), 'A')||strchr(subsDet[i][l].c_str(), 'T')||strchr(subsDet[i][l].c_str(), 'G')||
                //             strchr(subsDet[i][l].c_str(), 'C')||strchr(subsDet[i][l].c_str(), 'N')) ++e[j].trueMapCntSubs;
                //             else if (i==lineIdenFirst[n]&&l==lineIdenSec[n]){
                //                 ++e[j].misMapCntSubs;
                //                 ++n;
                //             }                     
                //             else if (strchr(softClip[i][l].c_str(), 'S')) ++e[j].softClipCntSubs;
                //             else if (softClip[i][l]=="*") ++e[j].noMapCntSubs;
                //            // else ++e[j].trueMapCntSubs;  
                //             ++m;
                //         }
                //         else{
                //             if (strchr(subsDet[i][l].c_str(), 'A')||strchr(subsDet[i][l].c_str(), 'T')||strchr(subsDet[i][l].c_str(), 'G')||
                //             strchr(subsDet[i][l].c_str(), 'C')||strchr(subsDet[i][l].c_str(), 'N')) {}//++e[j].trueMapCnt;
                //             else if (i==lineIdenFirst[n]&&l==lineIdenSec[n]){
                //                 ++e[j].misMapCnt;
                //                 ++n;
                //             }                     
                //             else if (strchr(softClip[i][l].c_str(), 'S')) ++e[j].softClipCnt;
                //             else if (softClip[i][l]=="*") ++e[j].noMapCnt;
                //             else ++e[j].trueMapCnt;       
                //         }
                //     }
                // }

                // if (readLen[i][l].length()==22){
                //                             // cout<<i<<" "<<l<<endl;
                //      if(readHead[i][l]==id[j]){
                //         ++f[j].totalReadCnt;
                //         if(l==subsLine[i][m]){
                //             if (strchr(subsDet[i][l].c_str(), 'A')||strchr(subsDet[i][l].c_str(), 'T')||strchr(subsDet[i][l].c_str(), 'G')||
                //             strchr(subsDet[i][l].c_str(), 'C')||strchr(subsDet[i][l].c_str(), 'N')) ++f[j].trueMapCntSubs;
                //             else if (i==lineIdenFirst[n]&&l==lineIdenSec[n]){
                //                 ++f[j].misMapCntSubs;
                //                 ++n;
                //             }                     
                //             else if (strchr(softClip[i][l].c_str(), 'S')) ++f[j].softClipCntSubs;
                //             else if (softClip[i][l]=="*") ++f[j].noMapCntSubs;
                //           //  else ++f[j].trueMapCntSubs;  
                //             ++m;
                //         }
                //         else{
                //             if (strchr(subsDet[i][l].c_str(), 'A')||strchr(subsDet[i][l].c_str(), 'T')||strchr(subsDet[i][l].c_str(), 'G')||
                //             strchr(subsDet[i][l].c_str(), 'C')||strchr(subsDet[i][l].c_str(), 'N')) {}//++f[j].trueMapCnt;
                //             else if (i==lineIdenFirst[n]&&l==lineIdenSec[n]){
                //                 ++f[j].misMapCnt;
                //                 ++n;
                //             }                     
                //             else if (strchr(softClip[i][l].c_str(), 'S')) ++f[j].softClipCnt;
                //             else if (softClip[i][l]=="*") ++f[j].noMapCnt;
                //             else ++f[j].trueMapCnt;       
                //         }
                //     }
                // }

            
        //     }
        //     ++j;
        // }  

        // cout<<a[1].totalReadCnt<<" "<<b[1].totalReadCnt<<" "<<c[1].totalReadCnt<<" "<<d[1].totalReadCnt<<" "<<e[1].totalReadCnt<<" "<<f[1].totalReadCnt;
        // for (i=0; i<totalRef; i++){
        //     for (j=0; j<30; j++){
        //         cout<<"hehe";
        //         cout<<subsLine[i][j]<<" ";
        //     }
        //     cout<<endl;
        // }
        // for (i=0; i<totalGff; i++){
        //     cout<<gffPosStrt[i]<<"\t"<<gffPosEnd[i]<<"\t"<<gffID[i]<<endl;
        // }

        // for (i=0; i<totalRef; i++){
        //     cout<<header[i]<<endl;
        //     // cout<<;
        //     cout<<"\t\t17\t18\t19\t20\t21\t22\t"<<endl;
        //     cout<<"TotalRead\t";
        //     cout<<a[i].totalReadCnt<<"\t"<<b[i].totalReadCnt<<"\t"<<c[i].totalReadCnt<<"\t"<<d[i].totalReadCnt<<"\t"<<e[i].totalReadCnt<<"\t"<<f[i].totalReadCnt<<endl;
        //     cout<<"NoSubsRead\t";
        //     cout<<a[i].noSubs<<"\t"<<b[i].noSubs<<"\t"<<c[i].noSubs<<"\t"<<d[i].noSubs<<"\t"<<e[i].noSubs<<"\t"<<f[i].noSubs<<endl;
        //     cout<<"softClip\t";
        //     cout<<a[i].softClipCnt<<"\t"<<b[i].softClipCnt<<"\t"<<c[i].softClipCnt<<"\t"<<d[i].softClipCnt<<"\t"<<e[i].softClipCnt<<"\t"<<f[i].softClipCnt<<endl;
        //     cout<<"correctMap\t";
        //     cout<<a[i].trueMapCnt<<"\t"<<b[i].trueMapCnt<<"\t"<<c[i].trueMapCnt<<"\t"<<d[i].trueMapCnt<<"\t"<<e[i].trueMapCnt<<"\t"<<f[i].trueMapCnt<<endl;
        //     cout<<"misMap\t\t";
        //     cout<<a[i].misMapCnt<<"\t"<<b[i].misMapCnt<<"\t"<<c[i].misMapCnt<<"\t"<<d[i].misMapCnt<<"\t"<<e[i].misMapCnt<<"\t"<<f[i].misMapCnt<<endl;
        //     cout<<"notMap\t\t";
        //     cout<<a[i].noMapCnt<<"\t"<<b[i].noMapCnt<<"\t"<<c[i].noMapCnt<<"\t"<<d[i].noMapCnt<<"\t"<<e[i].noMapCnt<<"\t"<<f[i].noMapCnt<<endl;

        //     cout<<"SubsRead\t";
        //     cout<<a[i].subs<<"\t"<<b[i].subs<<"\t"<<c[i].subs<<"\t"<<d[i].subs<<"\t"<<e[i].subs<<"\t"<<f[i].subs<<endl;
        //     cout<<"softClip\t";
        //     cout<<a[i].softClipCntSubs<<"\t"<<b[i].softClipCntSubs<<"\t"<<c[i].softClipCntSubs<<"\t"<<d[i].softClipCntSubs<<"\t"<<e[i].softClipCntSubs<<"\t"<<f[i].softClipCntSubs<<endl;
        //     cout<<"correctMap\t";
        //     cout<<a[i].trueMapCntSubs<<"\t"<<b[i].trueMapCntSubs<<"\t"<<c[i].trueMapCntSubs<<"\t"<<d[i].trueMapCntSubs<<"\t"<<e[i].trueMapCntSubs<<"\t"<<f[i].trueMapCntSubs<<endl;
        //     cout<<"misMap\t\t";
        //     cout<<a[i].misMapCntSubs<<"\t"<<b[i].misMapCntSubs<<"\t"<<c[i].misMapCntSubs<<"\t"<<d[i].misMapCntSubs<<"\t"<<e[i].misMapCntSubs<<"\t"<<f[i].misMapCntSubs<<endl;
        //     cout<<"notMap\t\t";
        //     cout<<a[i].noMapCntSubs<<"\t"<<b[i].noMapCntSubs<<"\t"<<c[i].noMapCntSubs<<"\t"<<d[i].noMapCntSubs<<"\t"<<e[i].noMapCntSubs<<"\t"<<f[i].noMapCntSubs<<endl<<endl;

        //     // cout<<a[i].softClipCnt<<"\t"<<b[i].softClipCnt<<"\t"<<c[i].softClipCnt<<"\t"<<d[i].softClipCnt<<"\t"<<e[i].softClipCnt<<"\t"<<f[i].softClipCnt<<endl;
        //     // // cout<<a[i].noSubsCnt<<"\t"<<b[i].noSubsCnt<<"\t"<<c[i].noSubsCnt<<"\t"<<d[i].noSubsCnt<<"\t"<<e[i].noSubsCnt<<"\t"<<f[i].noSubsCnt<<endl;
        //     // cout<<"SubsRead\t"<<a[i].subs<<"\t"<<b[i].subs<<"\t"<<c[i].subs<<"\t"<<d[i].subs<<"\t"<<e[i].subs<<"\t"<<f[i].subs<<"\t";
        //     // cout<<a[i].subsCnt+a[i].softClipCnt<<"\t"<<b[i].subsCnt+b[i].softClipCnt<<"\t"<<c[i].subsCnt+c[i].softClipCnt<<"\t"<<d[i].subsCnt+d[i].softClipCnt<<"\t"<<e[i].subsCnt+e[i].softClipCnt<<"\t"<<f[i].subsCnt+f[i].softClipCnt<<endl;
            
        //     // cout<<a[i].noMapCnt<<"\t"<<b[i].noMapCnt<<"\t"<<c[i].noMapCnt<<"\t"<<d[i].noMapCnt<<"\t"<<e[i].noMapCnt<<"\t"<<f[i].noMapCnt<<endl;   
        //     // cout<<"MisMap\t\t*\t*\t*\t*\t*\t*\t";        
        //     // cout<<a[i].misMapCnt<<"\t"<<b[i].misMapCnt<<"\t"<<c[i].misMapCnt<<"\t"<<d[i].misMapCnt<<"\t"<<e[i].misMapCnt<<"\t"<<f[i].misMapCnt<<endl<<endl;
        // }  

    // for (i=0; i<totalRef; i++){
    //     cout<<header[i]<<endl<<"--------------------"<<endl;
    //     cout<<"Total Read:\t"<<totalRead[i]<<"\t|Total Read:\t\t"<<totalReadCnt[i]<<"\t|"<<endl;
    //     cout<<"No Subs Read:\t"<<noSubsMap[i]<<"\t|No Subs Read Mapped:\t"<<noSubsCnt[i]<<"\t|"<<endl;
    //     cout<<"Subs Read:\t"<<subsMap[i]<<"\t|SoftClip:\t\t"<<softClipCnt[i]<<"\t|"<<endl;
    //     cout<<"\t\t\t|Subs Read Mapped:\t"<<softClipCnt[i]+subsCnt[i]<<"\t|"<<endl;
    //     cout<<"Num Subs:\t"<<subsBaseNum[i]<<"\t|read not mapped:\t"<<noMapCnt[i]<<"\t|"<<endl;
    //     cout<<"\t\t\t|incorrect mapped:\t"<<misMapCnt[i]<<"\t|"<<endl<<endl;
    // }

        // for (i=0; i< totalRef; i++){
        // cout<<a[i].totalReadRng<<"\t"<<b[i].totalReadRng<<"\t"<<c[i].totalReadRng<<"\t"<<d[i].totalReadRng<<"\t"<<e[i].totalReadRng<<"\t"<<f[i].totalReadRng<<endl;
        // cout<<a[i].noSubs<<"\t"<<b[i].noSubs<<"\t"<<c[i].noSubs<<"\t"<<d[i].noSubs<<"\t"<<e[i].noSubs<<"\t"<<f[i].noSubs<<endl;
        // cout<<a[i].subs<<"\t"<<b[i].subs<<"\t"<<c[i].subs<<"\t"<<d[i].subs<<"\t"<<e[i].subs<<"\t"<<f[i].subs<<endl<<endl;
        // }
//         return 0;
// }