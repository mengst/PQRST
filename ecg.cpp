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
	class Signal signal;
	signal.ReadFile("200.dat");
	double* data = signal.GetData(0);
	double* data_lead2 = signal.GetData(1);
	int size = signal.GetLength();
    double sr = signal.GetSR();
	class EcgAnnotation ann;
	
	
	vector<double> RRs;
	vector<double> QRS;
	int** qrsAnn = ann.GetQRS(data, size, sr, L"filters");         //get QRS complexes
	int** QrsAnnotation =ann.GetQrsAnnotation();
	cout << " GetQrsAnnotation(): " <<**QrsAnnotation  <<endl;
	wprintf(L" %d beats.\n", ann.GetQrsNumber());
	ann.GetEctopics(qrsAnn, ann.GetQrsNumber(), sr, &RRs, &QRS);        //label Ectopic beats
	int** ANN = ann.GetPTU(data, size, sr, L"filters", qrsAnn, ann.GetQrsNumber());     //find P,T waves
	cout << "ANN: " <<**ANN  <<endl;
	int annNum = 0;
	annNum = ann.GetEcgAnnotationSize();

	int startQRSms;
	int slutQRSms;
	int counter = 0;
	int counter_abnormal = 0;
	int R_msec = 0;//***********************
	int beat_nr = 0;//***********************
	int beat_nr_tri = 0;//**********************
	int beat1 = 0;//*******************
	int beat2 = 0;
	int beat3 = 0;
	int counter_ECT = 0;//******************
	int counter_ECT_t = 0;//*********************
	int counter_bigen=0;//************************
	int counter_bign_serie = 0;//***********************
	int counter_ves_serie = 0;//**********************
	int counter_rr = 0;//*******************************
	int buf_msec = 0;//***************************
	int buf_msec_dif = 0;//*************************
	//int buf_beat_ves = 0;							
	int flag_abnormal = 0;//************************
	int flag_r = 0; int flag_p = 0; int flag_t = 0; int flag_n = 0; int flag_e = 0;//*************
	int flag_ves = 0; int flag_no_p = 0; int flagga_ves_s = 0;//************************
	int flagga_big = 0; int flag_r_str = 0; int beat_r_str = 0;//*************************
	int flag_r_stp = 0; int beat_r_stp = 0; int r_sta_ms = 0; int r_stp_ms = 0;
	float svt = 0;
	int cntr = 0;//??????????????????????????
	int cntr_u = 0;
	int sum_rr = 0;//****************************
	int rr = 0;//*********************
	int medel_rr = 0;//**********************
	wchar_t startECT = 'E';//****************************
	wchar_t R_comp = 'R';//******************************
	wchar_t P_wave = 'p';//********************************
	wchar_t T_wave = 'T';//***********************************
	wchar_t startA = 'A';
	wchar_t startQRS = 'N';
	wchar_t r_sta = 'L';
	wchar_t r_stp = 'K';
	wchar_t slutQRS = ')';
	wchar_t U = 'U';
	wchar_t r = 'r';
	vector<double> BBs;
	double QRs;
	int h, m, s, ms;
	int msec = int(((double)size / sr) * 1000.0);

	for (int i = 0; i < annNum; i++) { 
									

	int smpl = ANN[i][0];
	int type = ANN[i][1];
									

	msec = int(((double)smpl / sr) * 1000.0);
	signal.mSecToTime(msec, h, m, s, ms);
	wchar_t * test = anncodes[type];
							
									
				
//		if (i < 1000){ wprintf(L"   %04d %02d:%02d:%02d %03d  %s    \n", smpl, h, m, s, ms, anncodes[type]); }

	if (*test == r_sta)// when rytm changes to shorter set the flag to 1 for counting the beats
	{
		flag_r_str = 1;
	}

	if (*test == r_stp && flag_r_str==1)// the rytm changes to longer
	{
		flag_r_stp = 1;
	}
							


	if (*test == P_wave)// when a P-wave detects if a beat, the p-flag will set to 1
	{
		flag_p = 1;
	}

	if ((*test == r_sta) || (*test == r_stp)||(*test == startECT) || (*test == startQRS) || (*test == startA) || (*test == U))// sparar beat nr när den hittar en R_våg
	{	
										
	//	if (beat_nr>1 && beat_nr<200/*ann.GetQrsNumber()-2*/){ wprintf(L"beatnr:%d Pvåg:%d  %04d %02d:%02d:%02d.%05d  %s  rrs:%.2f  QRS:%.3f  \n", beat_nr, flag_p, smpl, h, m, s, ms, anncodes[type],RRs[beat_nr-1],QRS[beat_nr]); }
		beat_nr++;
										

		if (flag_r_str == 1)// analyze if a rytm change has begun and saves the time
		{ 
			if (beat_r_str == 0){ r_sta_ms = msec; }
			beat_r_str++;// counts the number of beats when rythm change occures
									
		}
										
		// counts the beats between rythm changes and determines if VT
		if (flag_r_str == 1 && flag_r_stp == 1)
		{ 
											
			r_stp_ms = int(((double)ANN[i-1][0] / sr) * 1000.0);
			svt =(r_stp_ms - r_sta_ms);
			svt=(((beat_r_str-1)*60) / svt)*1000;
			if (beat_r_str>3 && svt > 100){
				wprintf(L" \n SVT över 100\n");
				flag_abnormal = 1;
			}

			flag_r_str = 0;
			flag_r_stp = 0;
			beat_r_str = 0;
			svt = 0;
		}
	}

	if (*test == startQRS)
	{
		flag_n = 1;
	}

									

	//****************hittar slag utan p våg****************************
	if (((*test == startQRS) && (flag_p == 0))/* || (*test == r_stp && flag_p == 0)*/)
	{
	//	wprintf(L"no p msec %d %10d %02d:%02d:%02d.%03d   %s\n", msec, smpl, h, m, s, ms, anncodes[type]);
		if (cntr_u > 1) { flag_no_p = 1; }
		cntr_u++;
										
	}

									

	//*********** hittar extra slag, bigenimi, i serie, avgör patologiskt*****
	if ((*test == startECT) || (*test == U) || (*test == startA) || flag_no_p == 1 || (*test == r_sta)){
		counter_ECT++;
		counter_ECT_t++;
		//wprintf(L"no p msec %d %10d %02d:%02d:%02d.%03d   %s\n", msec, smpl, h, m, s, ms, anncodes[type]);
		flag_e = 1;
		flag_no_p = 0;

		if (counter_ECT_t == 1)// when first ect beat is found it will search for more in two minutes
		{
			buf_msec = msec + 120000;
			buf_msec_dif = msec;
		}


		if ((msec <= buf_msec) && (counter_ECT_t > 5))// om det förekommer fler än 5 slagpå 2min
		{
		//	wprintf(L"abnormal/ många extraslag på 2 min\n");
		//	wprintf(L"%d %10d %02d:%02d:%02d.%03d   %s\n", msec, smpl, h, m, s, ms, anncodes[type]);

			flag_abnormal = 1;
			counter_ECT_t = 0;
			counter_abnormal++;

		}

		if ((counter_ECT_t <= 5) && (msec >= buf_msec))
		{
			counter_ECT_t = 0;
		//	wprintf(L"få extra slag under fem min, nollställs\n");
		}


		// bigenemi detects start of bigenemi
										
		if ( beat1 ==(beat2 + 2) && beat_nr!=(beat1+1) && beat2!=(beat3+1)&& flagga_ves_s==0)// sets flag big to 1
		{
		//	wprintf(L" beatNR: %d  beat1: %d beat2:%d  beat3:%d \n", beat_nr, beat1, beat2, beat3);
			flagga_big = 1;
		//	wprintf(L"big slag %d %10d %02d:%02d:%02d.%03d   %s \n\n", msec, smpl, h, m, s, ms, anncodes[type]);
											
											
		}

		// counts nr of beginemi and finds when bigenemi stops
										
		if ((beat_nr !=beat1) && (flagga_big == 1) && (beat1 != (beat2 + 2)))// counts big i serie
		{
		//	wprintf(L" big slut: %d  beat1: %d beat2:%d  beat3:%d \n", beat_nr, beat1, beat2, beat3);
			counter_bigen++;
		//	wprintf(L"bigenemi %d \n\n", counter_bigen);
			flagga_big = 0;
		//	wprintf(L"msec avvikande slag %d %10d %02d:%02d:%02d.%03d   %s\n", msec, smpl, h, m, s, ms, anncodes[type]);

		}

        // VES i serie *************************
		if (beat_nr == (beat1 + 1))
		{	
			flagga_ves_s = 1;
	//		wprintf(L"msec avvikande slag %d %10d %02d:%02d:%02d.%03d   %s\n", msec, smpl, h, m, s, ms, anncodes[type]);
											
	//		wprintf(L"ves flag %d \n", flagga_ves_s);
			//if (beat1 = beat2 + 1) { counter_bigen--; }
		}
	// slut trigenemi****************

		if (beat_nr != (beat1 +1) && flagga_ves_s == 1)// sets flag big to 1
		{
			counter_ves_serie++;
			flagga_ves_s=0;
			flag_ves = 1;

		}

		// check if a singel beat is after ves i serie
										

	/*	if (beat_nr==(buffer_beat_nr+1))
		{
			flagga_big = 1;
			counter_ves_serie = 0;// flagga för trigenemi sätts till noll
			if ((flagga_big==1) && (beat_nr_tri == (beat_nr - 2)))// trigenemi
			{
			//	wprintf(L"trigenemi i  \n");
			//	wprintf(L"%d %10d %02d:%02d:%02d.%03d   %s\n", msec, smpl, h, m, s, ms, anncodes[type]);
				counter_ves_serie = 1;
				flag_abnormal = 1;
				counter_abnormal++;
				flagga_big = 0;
			}

			beat_nr_tri = buffer_beat_nr;// sparar beatnr för att avgöra trigenemi
		//	wprintf(L"bigenemi\n");
		//	wprintf(L"%d %10d %02d:%02d:%02d.%03d   %s\n",msec, smpl, h, m, s, ms, anncodes[type]);
			counter_bigen++;


		}*/

		if ((counter_bigen > 2) && (counter_ves_serie>1))// detekterar bigenemi i serie
		{
			counter_bign_serie++;
				flag_abnormal = 1;
				counter_abnormal++;
		//		wprintf(L"bigenemi i serie %d \n", counter_bign_serie);
		//		wprintf(L"%d %10d %02d:%02d:%02d.%03d   %s\n", msec, smpl, h, m, s, ms, anncodes[type]);

		}


		startQRSms = msec;
		if (beat_nr != beat1){
			beat3 = beat2;
			beat2 = beat1;
			beat1 = beat_nr;
		}
										
										
	/*	buf_beat_ves = buffer_beat_nr;
		buffer_beat_nr = beat_nr;*/
	}

	//************************slut på ECT funktion******************************


									



	if (*test == startQRS){

		//flag_n = 1;
		//	wprintf(L"msec %d %10d %02d:%02d:%02d.%03d   %s\n", msec, smpl, h, m, s, ms, anncodes[type]);
		startQRSms = msec;
	}

/*	if (*test == T_wave)
	{
			flag_t = 1;

		}*/

	if (*test == R_comp)
	{
		flag_r = 1;

	}

		if (*test == slutQRS){
		//	wprintf(L"msec %d %10d %02d:%02d:%02d.%03d   %s\n", msec, smpl, h, m, s, ms, anncodes[type]);
			slutQRSms = msec;
			//wprintf(L"QRS time: %d\n",(slutQRSms-startQRSms));
		}
		if((slutQRSms-startQRSms)> 120){
			//counter++;
		//	wprintf(L"QRS time: %d\n",(slutQRSms-startQRSms));
		}


		//************************ flaggor för beats nollställs när ett r våg upptäcks************
		if (flag_r == 1)
		{
			flag_n = 0;
			flag_p = 0;
			flag_e = 0;
			flag_r = 0;
			flag_t = 0;
			flag_no_p = 0;
		}

}

}
/*void main(){
	Magnus
	class Signal signal;
	signal.ReadFile("200.dat");
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
	vector<int> qtInterval;
	int pWave = 0;
	int tWave = 0;
	int qt = 0;
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
		while(ANN[j+pk][1] != 45 && ((j+pk) < ann.GetEcgAnnotationSize()-1) && pk < 9){	//Find the end of the T wave, which is equal to 45. (might add pk<5) 
				pk++;
		}
		if(ANN[j+pk][1] == 45){		//if QT interval is found
			qt = ANN[(j+pk)][0]-ANN[j][0];	//calculate the length of the wave, the result is in samples
		}
		else{
			qt = 0;	//if no T-wave is found
		}
		qtInterval.push_back(qt);
		pk = 0;

		
		while(ANN[j+pk][1] != 44 && ((j+pk) < ann.GetEcgAnnotationSize()-1) && pk < 7){	//Find the start of the T wave, which is equal to 44. (might add pk<5) 
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
	int pWaveMilliSec, tWaveMilliSec, qtInt;
	
	//Loop to print the information for each heart beat
	for (int i = 0; i < ann.GetQrsNumber(); i++) {

		int smpl = qrsAnn[2*i][0];
		int type = qrsAnn[(2*i)][1];
		int qrsComplex = ((qrsAnn[(2*i)+1][0])-(qrsAnn[2*i][0]));

		//convert to milliseconds
        qrsComplex = int(((double)qrsComplex / sr) * 1000.0);
		pWaveMilliSec = int(((double)pWaveLength[i] / sr) * 1000.0);
		tWaveMilliSec = int(((double)tWaveLength[i] / sr) * 1000.0);
		qtInt = int(((double)qtInterval[i] / sr) * 1000.0);
		msec = int(((double)smpl / sr) * 1000.0);
        signal.mSecToTime(msec, h, m, s, ms);
		
		//just the first heart beat
		if(i==0){
			wprintf(L"%10d %02d:%02d:%02d.%03d RR.int: 0.000 T wave: %d QT: %d %s\n", smpl, h, m, s, ms,tWaveMilliSec, qtInt, anncodes[type] );
		}
		
		//If a normal heart beat, but the QRS complex exceeds 120 ms
		//if (type == 1 && qrsComplex >= 126 && (i-1)>=0){
		//	wprintf(L"%10d %02d:%02d:%02d.%03d RR.int: %.3f QRS complex: %d PACE\n", smpl, h, m, s, ms, RRs[i-1], qrsComplex );
		//}

		//Print all other heart beats
		if((i-1)>=0 ){
			wprintf(L"%10d %02d:%02d:%02d.%03d RR.int: %.3f P wave: %d T wave: %d QT: %d QRS: %d %s\n\n", smpl, h, m, s, ms, RRs[i-1], pWaveMilliSec, tWaveMilliSec, qtInt, qrsComplex, anncodes[type] );
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
	
	 
}*/
	

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