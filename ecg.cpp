/*

ecg.cpp - ECG annotation console app based lib.cpp ECG library.
Copyright (C) 2007 YURIY V. CHESNOKOV

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


You may contact the author by e-mail (chesnokov.yuriy@gmail.com or chesnokov_yuriy@mail.ru)
or postal mail (Unilever Centre for Molecular Sciences Informatics,
University Chemical Laboratory, Cambridge University, Lensfield Road, Cambridge, CB2 1EW, UK)

*/


#include "stdafx.h"
#include "lib\lib.h"

void tic();
void toc();
void help();
void change_extension(wchar_t* path, const wchar_t* ext);
int parse_params(class EcgAnnotation &ann);

wchar_t params[_MAX_PATH] = L"params";

void main(){
	const wchar_t * argv;
	argv = L"106.dat";
	class Signal signal;
	signal.ReadFile("106.dat");
	double* data = signal.GetData(0);
	double* data_lead2 = signal.GetData(1);
	int size = signal.GetLength();
    double sr = signal.GetSR();
	class EcgAnnotation ann;
	int** qrsAnn = ann.GetQRS(data, size, sr, L"filters");         //get QRS complexes
	
	vector<double> RRs;
	ann.GetEctopics(qrsAnn, ann.GetQrsNumber(), sr, &RRs);        //label Ectopic beats', RRs = RR avstånd
	int h, m, s, ms;
    int msec = int(((double)size / sr) * 1000.0);
	int annNum = 0;
	annNum = ann.GetEcgAnnotationSize();
	
	
	int tt = 0;
	wprintf(L"Total beats: %d\n", ann.GetQrsNumber()); 
	
	
	
    int** ANN = ann.GetPTU(data, size, sr, L"filters", qrsAnn, ann.GetQrsNumber());     //find P,T waves
	
	vector<int> pWaveLength;
	vector<int> tWaveLength;
	int pWave = 0;
	int tWave = 0;
	int j =0;
	int smpl = 0;
	int pk =0;
	pWaveLength.push_back(0);		//First beat doesn't have a P-wave

	//Loop to find P- and T-wave for each beat
	for (int i = 0; i < ann.GetQrsNumber(); i++){
		smpl = qrsAnn[2*i][0];	//The sample for the start of every QRS Complex (N and ECG)		
		
		while(smpl != ANN[j][0]){	//Find the place for this sample in the array where all annotations are stored (ANN)
			j++;
		}
		while(ANN[j+pk][1] != 44 && (j+pk) < ann.GetEcgAnnotationSize()-1){	//Find the start of the T wave, which is equal to 44. (might add pk<5) 
			pk++;
		}
		if(ANN[j+pk][1] == 44){		//if T-wave is found
			tWave = ANN[(j+pk)+2][0]-ANN[j+pk][0];	//calculate the length of the wave, the result is in samples
		}
		else{
			tWave = 0;	//if no T-wave is found
		}
		tWaveLength.push_back(tWave);	//Add the P-wave length to an array
		pk =0;

		//Find the P-wave
		if(j>3){	
			if(ANN[j-3][1] == 42){	//If the 3rd place before the start of the QRS complex is indicating a start of a p-wave  
				pWave = ANN[j-1][0] - ANN[j-3][0];	//calculate the length of the wave, the result is in samples
			}
			else{
				pWave = 0;
			}
			pWaveLength.push_back(pWave);
		}
		
	}
	int pWaveMilliSec, tWaveMilliSec;
	
	//Loop to print the information for each heart beat
	for (int i = 0; i < ann.GetQrsNumber(); i++) {

		int smpl = qrsAnn[2*i][0];
		int type = qrsAnn[(2*i)][1];
		int qrsComplex = ((qrsAnn[(2*i)+1][0])-(qrsAnn[2*i][0]));

		//convert to milliseconds
        qrsComplex = int(((double)qrsComplex / sr) * 1000.0);
		pWaveMilliSec = int(((double)pWaveLength[i] / sr) * 1000.0);
		tWaveMilliSec = int(((double)tWaveLength[i] / sr) * 1000.0);
		msec = int(((double)smpl / sr) * 1000.0);
        signal.mSecToTime(msec, h, m, s, ms);
		
		//just the first heart beat
		if(i==0){
			wprintf(L"%10d %02d:%02d:%02d.%03d RR.int: 0.000 %s\n", smpl, h, m, s, ms,anncodes[type] );
		}
		
		//If a normal heart beat, but the QRS complex exceeds 120 ms
		if (type == 1 && qrsComplex >= 120 && (i-1)>=0){
			wprintf(L"%10d %02d:%02d:%02d.%03d RR.int: %.3f QRS complex: %d PACE\n", smpl, h, m, s, ms, RRs[i-1], qrsComplex );
		}

		//Print all other heart beats
		if((i-1)>=0 && (type == 1 || type == 46) ){
			wprintf(L"%10d %02d:%02d:%02d.%03d RR.int: %.3f P wave: %d T wave: %d  QRS: %d %s\n", smpl, h, m, s, ms, RRs[i-1], pWaveMilliSec, tWaveMilliSec, qrsComplex, anncodes[type] );
		}
    }

	annNum = ann.GetEcgAnnotationSize();	//get the number of all annotations

	//print information about all annotations
	for (int i = 0; i < annNum; i++) {
		int sampl = ANN[i][0];
        int type = ANN[i][1];

        msec = int(((double)sampl / sr) * 1000.0);
        signal.mSecToTime(msec, h, m, s, ms);

        wprintf(L"%10d %02d:%02d:%02d.%03d   %s\n", sampl, h, m, s, ms, anncodes[type]);
	}
	
	 
}
	

void help()
{
        wprintf(L"usage: ecg.exe physionetfile.dat [LeadNumber] [params]\n");
        wprintf(L"       do not forget about \\filters dir to be present.");
}

static LARGE_INTEGER m_nFreq;
static LARGE_INTEGER m_nBeginTime;

void tic()
{
        QueryPerformanceFrequency(&m_nFreq);
        QueryPerformanceCounter(&m_nBeginTime);
}
void toc()
{
        LARGE_INTEGER nEndTime;
        __int64 nCalcTime;

        QueryPerformanceCounter(&nEndTime);
        nCalcTime = (nEndTime.QuadPart - m_nBeginTime.QuadPart) * 1000 / m_nFreq.QuadPart;

        wprintf(L" processing time: %d ms\n", nCalcTime);
}

void change_extension(wchar_t* path, const wchar_t* ext)
{
        for (int i = (int)wcslen(path) - 1; i > 0; i--) {
                if (path[i] == '.') {
                        path[i] = 0;
                        wcscat(path, ext);
                        return;
                }
        }
        wcscat(path, ext);
}

int parse_params(class EcgAnnotation &ann)
{
        FILE* fp = _wfopen(params, L"rt");
        if (fp != 0) {
                ANNHDR hdr;
                int res = 0;
                res = fwscanf(fp, L"%*s %d %*s %d"
                                  L"%*s %lf %*s %lf %*s %lf %*s %d %*s %lf"
                                  L"%*s %lf %*s %lf %*s %lf %*s %lf"
                                  L"%*s %lf %*s %lf %*s %d",  
                                   &hdr.minbpm, &hdr.maxbpm,
                                   &hdr.minQRS, &hdr.maxQRS, &hdr.qrsFreq, &hdr.ampQRS, &hdr.minUmV,
                                   &hdr.minPQ, &hdr.maxPQ, &hdr.minQT, &hdr.maxQT, 
                                   &hdr.pFreq, &hdr.tFreq, &hdr.biTwave);
                if (res == 14) {
                        PANNHDR phdr = ann.GetAnnotationHeader();
                        memcpy(phdr, &hdr, sizeof(ANNHDR));
                        wprintf(L" using annotation params from file %s\n", params);
                        wprintf(L"  minBpm  %d\n"  
                                L"  maxBpm  %d\n" 
                                L"  minQRS  %lg\n" 
                                L"  maxQRS  %lg\n" 
                                L" qrsFreq  %lg\n"
                                L"  ampQRS  %d\n"
                                L"  minUmV  %lg\n"
                                L"   minPQ  %lg\n"
                                L"   maxPQ  %lg\n" 
                                L"   minQT  %1f\n" 
                                L"   maxQT  %lg\n" 
                                L"   pFreq  %lg\n"  
                                L"   tFreq  %lg\n"
                                L" biTwave  %d\n\n", hdr.minbpm, hdr.maxbpm,
                                                   hdr.minQRS, hdr.maxQRS, hdr.qrsFreq, hdr.ampQRS, hdr.minUmV,
                                                   hdr.minPQ, hdr.maxPQ, hdr.minQT, hdr.maxQT, 
                                                   hdr.pFreq, hdr.tFreq, hdr.biTwave);
                        fclose(fp);
                        return 0;
                }
                else {
                        fclose(fp);
                        wprintf(L" failed to read %s annotation params file, using default ones instead.\n", params);
                        return res;                                
                }
        }
        else {
                wprintf(L" failed to open %s annotation params file, using default ones instead.\n", params);
                return -1;
        }
}