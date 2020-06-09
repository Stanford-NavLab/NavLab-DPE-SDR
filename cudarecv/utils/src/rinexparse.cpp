
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <algorithm> // std::find
#include "rinexparse.h"
#include "auxil.h"
#include "errorhandler.h"

/**
 *  Adapted from rtklib.h function readrnx and functions called within
 *  Expected that all arguments are initialized prior to function call
 *
 *  args	:	FILE *fp				Pointer to RINEX file
 *  			navlite_t *nav			Pointer to array of navs
 *  			g_time target_eph_time	Target time for ephemerides
 *
 */
int
dsp::utils::ReadRinex(FILE *fp, const char *opt, std::vector<dsp::utils::ephSet_t> &nav) {

	char filebuf[1024];

	char type=' ';
	double ver = 2.10;
	int sys,tsys=dsp::utils::TSYS_GPS;
	dsp::utils::rinexHeaderParams_t* curParams = new dsp::utils::rinexHeaderParams_t(); // This doesn't need to live outside one RINEX file

	dsp::utils::gtime_t t = {0}; // Observation time (not sure if this is useful in any way...)


	// Read RINEX header
	if (ReadRinexHeader(fp, &ver, &type, &sys, &tsys, nav, curParams) != 0) {
		funcLogHelperM("Failed to read RINEX header");
		delete[] curParams;
		return -1;
	}

	// Process RINEX body
    switch (type) {
        //case 'O': return readrnxobs(fp,ts,te,tint,opt,index,ver,&tsys,tobs,obs, sta);
        case 'N': return ReadRinexBody(fp, opt, ver, sys, nav, curParams);
        //case 'G': return readrnxnav(fp,opt,ver,SYS_GLO,nav);
        //case 'H': return readrnxnav(fp,opt,ver,SYS_SBS,nav);
        //case 'J': return readrnxnav(fp,opt,ver,SYS_QZS,nav); /* extension */
        //case 'L': return readrnxnav(fp,opt,ver,SYS_GAL,nav); /* extension */
        //case 'C': return readrnxclk(fp,opt,index,nav);
        default:
        	funcLogHelperM("Unsupported RINEX body type");
        	delete[] curParams;
        	return -1;
        	break;
	}


    delete[] curParams;
	return 0;
}
namespace dsp { namespace utils {
	namespace {
		/**
		 *  Adapted from readrnxfp in rinex.c in rtklib
		 */
		int ReadRinexHeader(FILE *fp, double *ver, char *type, int *sys, int *tsys, std::vector<dsp::utils::ephSet_t> &nav, dsp::utils::rinexHeaderParams_t *curParams) {
			double bias;
			char buff[utils::MAXRNXLEN],*label=buff+60;
			int i=0, block=0, sat;

			*ver=2.10, *type=' '; *sys=SYS_GPS;

			// Read in a buffer's worth (or until a newline/EOF is hit)
			while (fgets(buff,MAXRNXLEN,fp)) {

				// Less than 60 chars is an incomplete line in RINEX; nothing valuable here
				if (strlen(buff)<=60) continue;

				// Found the file type; set params
				else if (strstr(label,"RINEX VERSION / TYPE")) {
					*ver=auxil::str2num(buff,0,9);
					*type=*(buff+20);

					/* satellite system (code for commented cases not included here) */
					switch (*(buff+40)) {
						case ' ':
						case 'G': *sys=SYS_GPS;  *tsys=TSYS_GPS; break;
						//case 'R': *sys=SYS_GLO;  *tsys=TSYS_UTC; break;
						//case 'E': *sys=SYS_GAL;  *tsys=TSYS_GAL; break; /* v.2.12 */
						//case 'S': *sys=SYS_SBS;  *tsys=TSYS_GPS; break;
						//case 'J': *sys=SYS_QZS;  *tsys=TSYS_QZS; break; /* v.3.02 */
						//case 'C': *sys=SYS_CMP;  *tsys=TSYS_CMP; break; /* v.2.12 */
						//case 'I': *sys=SYS_IRN;  *tsys=TSYS_IRN; break; /* v.3.03 */
						//case 'M': *sys=SYS_NONE; *tsys=TSYS_GPS; break; /* mixed */
						default :
							funcLogHelperM("RINEX satellite system " << *(buff+40) << " unsupported");
							return -1;
					}
					continue;
				}
				else if (strstr(label,"PGM / RUN BY / DATE")) continue;
				else if (strstr(label,"COMMENT")) {
					// Include other if-elses from readrnxh to parse other params
					continue;
				}

				// Determine filetype (all prior if-else statements have continues,
				// so the program only gets here if all prior checks fail)
				switch (*type) {
					//case 'O': decode_obsh(fp,buff,*ver,tsys,tobs,nav,sta); break;
					case 'N': DecodeRinexHeader (buff,curParams); break;
					//case 'G': decode_gnavh(buff,nav); break;	/* currently unsupported */
					//case 'H': decode_hnavh(buff,nav); break;	/* currently unsupported */
					case 'J': DecodeRinexHeader (buff,curParams); break; 	/* extension */
					case 'L': DecodeRinexHeader (buff,curParams); break; 	/* extension */
					default:
						funcLogHelperM("RINEX filetype " << *type << " unsupported");
						return -1;
						break;
				}
				if (strstr(label,"END OF HEADER")) return 0;

				if (++i>=MAXPOSHEAD&&*type==' ') {
					funcLogHelperM("received non-RINEX file");
					return -1; /* no rinex file */
				}
			}


			return 0;
		}

		/**
		 *  Taken from decode_navh in rinex.c in rtklib
		 */
		int DecodeRinexHeader(char *buff, dsp::utils::rinexHeaderParams_t *curParams) {

			int i,j;
			char *label=buff+60;

			if      (strstr(label,"ION ALPHA"           )) { /* opt ver.2 */
				for (i=0,j=2;i<4;i++,j+=12) curParams->ion_gps[i]=auxil::str2num(buff,j,12);
			}
			else if (strstr(label,"ION BETA"            )) { /* opt ver.2 */
				for (i=0,j=2;i<4;i++,j+=12) curParams->ion_gps[i+4]=auxil::str2num(buff,j,12);
			}
			else if (strstr(label,"DELTA-UTC: A0,A1,T,W")) { /* opt ver.2 */
				for (i=0,j=3;i<2;i++,j+=19) curParams->utc_gps[i]=auxil::str2num(buff,j,19);
				for (;i<4;i++,j+=9) curParams->utc_gps[i]=auxil::str2num(buff,j,9);
			}
			else if (strstr(label,"IONOSPHERIC CORR"    )) { /* opt ver.3 */
				if (!strncmp(buff,"GPSA",4)) {
					for (i=0,j=5;i<4;i++,j+=12) curParams->ion_gps[i]=auxil::str2num(buff,j,12);
				}
				else if (!strncmp(buff,"GPSB",4)) {
					for (i=0,j=5;i<4;i++,j+=12) curParams->ion_gps[i+4]=auxil::str2num(buff,j,12);
				}
				// Include other if-elses from decode_navh to parse other RINEX files
			}
			else if (strstr(label,"TIME SYSTEM CORR"    )) { /* opt ver.3 */
				if (!strncmp(buff,"GPUT",4)) {
					curParams->utc_gps[0]=auxil::str2num(buff, 5,17);
					curParams->utc_gps[1]=auxil::str2num(buff,22,16);
					curParams->utc_gps[2]=auxil::str2num(buff,38, 7);
					curParams->utc_gps[3]=auxil::str2num(buff,45, 5);
				}
				// Include other if-elses from decode_navh to parse other RINEX files
			}
			else if (strstr(label,"LEAP SECONDS"        )) { /* opt */
				curParams->leaps=(int)auxil::str2num(buff,0,6);
			}

			return 0;
		}

		/**
		 *  Importing the RINEX body data is broken down into 3 functions:
		 *  	1. DecodeRinexBody 	fills an ephSet_t with appropriate data
		 *  	2. ParseRinexBody	breaks down strings read from the file into values for DecodeRinexBody
		 *  	3. ReadRinexBody	iterates over ParseRinexBody until the end of the file, placing ephSet_t's accordingly
		 *
		 */


		/**
		 *  Adapated from readrnxnav from rinex.c in rtklib
		 */
		int ReadRinexBody(FILE *fp, const char *opt, double ver, int sys, std::vector<dsp::utils::ephSet_t> &nav, dsp::utils::rinexHeaderParams_t *curParams) {

			dsp::utils::eph_t curEph;
			int stat, type;
			dsp::utils::ephSet_t dummyEph;
			int toeIdx;

			// Parse and store the RINEX data until the EOF
			while((stat=ParseRinexBody(fp, opt, ver, sys, &type, &curEph))>=0) {

				// To support more GNSS's, add if-else here to choose between the following code (add_eph adaptation)
				// and other new code (add_geph, add_seph adaptations)

				// Check if this toe is stored in the set of ephemerides, and make space for it if not
				// (Lambda black magic sorcery taken from here: http://www.cplusplus.com/forum/beginner/169178/
				int tempToes = curEph.toes;
				auto toeIt = std::find_if(
						nav.begin(),
						nav.end(),
						[tempToes](const dsp::utils::ephSet_t& nav)
							{return nav.ephToes == tempToes;}
				);
				toeIdx = toeIt - nav.begin();

				if (toeIdx == nav.size()) {
					nav.push_back(dsp::utils::ephSet_t());
					nav[toeIdx].ephToes = curEph.toes;
					std::copy(std::begin(nav[toeIdx].utc_gps), std::end(nav[toeIdx].utc_gps), std::begin(curParams->utc_gps));
					std::copy(std::begin(nav[toeIdx].ion_gps), std::end(nav[toeIdx].ion_gps), std::begin(curParams->ion_gps));
					nav[toeIdx].leaps 	= curParams->leaps;
				}

				// Store curEph accordingly, set corresponding valid flag
				nav[toeIdx].eph[curEph.sat-1] = curEph;
				nav[toeIdx].ephValid[curEph.sat-1] = 1;
			}

			return 0;
		}

		/**
		 *  Adapted from readrnxnavb from rinex.c in rtklib
		 */
		int ParseRinexBody(FILE *fp, const char *opt, double ver, int sys, int *type, dsp::utils::eph_t *eph) {

			dsp::utils::gtime_t toc;
			double data[64];
			int i=0,j,prn,sat=0,sp=3,mask;
			char buff[MAXRNXLEN],id[8]="",*p;

			mask=set_sysmask(opt);

			while (fgets(buff,MAXRNXLEN,fp)) {

				if (i==0) {

					/* decode satellite field */
					if (ver>=3.0||sys==SYS_GAL||sys==SYS_QZS) { /* ver.3 or GAL/QZS */
						//strncpy(id,buff,3);
						//sat=satid2no(id);
						//sp=4;
						//if (ver>=3.0) sys=satsys(sat,NULL);
						funcLogHelperM("Received unsupported RINEX version " << ver << ". Should be < 3.0.");
					}
					else {
						prn=(int)auxil::str2num(buff,0,2);

						if (sys==SYS_SBS) {
							sat=satno(SYS_SBS,prn+100);
						}
						else if (sys==SYS_GLO) {
							sat=satno(SYS_GLO,prn);
						}
						else if (93<=prn&&prn<=97) { /* extension */
							sat=satno(SYS_QZS,prn+100);
						}
						else sat=satno(SYS_GPS,prn);
					}
					/* decode toc field */
					if (dsp::utils::str2time(buff+sp,0,19,&toc)) {
						funcLogHelperM("rinex nav toc error: " << buff);
						return 0;
					}
					/* decode data fields */
					for (j=0,p=buff+sp+19;j<3;j++,p+=19) {
						data[i++]=auxil::str2num(p,0,19);
					}
				}
				else {
					/* decode data fields */
					for (j=0,p=buff+sp;j<4;j++,p+=19) {
						data[i++]=auxil::str2num(p,0,19);
					}
					/* decode ephemeris */
					if (sys==SYS_GLO&&i>=15) {
						if (!(mask&sys)) return 0;
						//*type=1;
						//return decode_geph(ver,sat,toc,data,geph);
						funcLogHelperM("unsupported ephemeris type");
					}
					else if (sys==SYS_SBS&&i>=15) {
						if (!(mask&sys)) return 0;
						//*type=2;
						//return decode_seph(ver,sat,toc,data,seph);
						funcLogHelperM("unsupported ephemeris type");
					}
					else if (i>=31) {
						if (!(mask&sys)) return 0;
						//*type=0;
						return DecodeRinexBody(ver,sat,toc,data,eph);
					}
				}
			}
			return -1;
		}


		/**
		 *  Adapted from decode_eph from rinex.c in rtklib
		 */
		int DecodeRinexBody(double ver, int sat, dsp::utils::gtime_t toc, const double *data, dsp::utils::eph_t *eph) {

			//dsp::utils::eph_t eph0={0}; // Don't need/use a local eph with the way ReadRinexBody is written
			int sys;

			sys=satsys(sat,NULL);

			//if (!(sys&(SYS_GPS|SYS_GAL|SYS_QZS|SYS_CMP|SYS_IRN))) {
			if (!(sys&(SYS_GPS))) { // Only GPS is currently supported; abort if other satellite system
				funcLogHelperM("ephemeris error: invalid satellite sat: " << sat);
				return 0;
			}
			//*eph=eph0; // Don't need/use a local eph with the way ReadRinexBody is written

			eph->sat=sat;
			eph->toc=toc;
			eph->tocs=dsp::utils::time2gpst(toc, &eph->week); // Added to avoid time conversion during kernel functions

			eph->f0=data[0];
			eph->f1=data[1];
			eph->f2=data[2];

			eph->sqrt_A=data[10]; eph->e=data[ 8]; eph->i0  =data[15]; eph->OMG0=data[13];
			eph->omg =data[17]; eph->M0 =data[ 6]; eph->deln=data[ 5]; eph->OMGd=data[18];
			eph->idot=data[19]; eph->crc=data[16]; eph->crs =data[ 4]; eph->cuc =data[ 7];
			eph->cus =data[ 9]; eph->cic=data[12]; eph->cis =data[14];

			eph->A=(eph->sqrt_A*eph->sqrt_A);
			eph->e_sqr = (eph->e * eph->e);

			if (sys==SYS_GPS||sys==SYS_QZS) {
				eph->iode=(int)data[ 3];      /* IODE */
				eph->iodc=(int)data[26];      /* IODC */
				eph->toes=     data[11];      /* toe (s) in gps week */
				eph->week=(int)data[21];      /* gps week */
				eph->toe=dsp::utils::adjweek(dsp::utils::gpst2time(eph->week,data[11]),toc);
				eph->ttr=dsp::utils::adjweek(dsp::utils::gpst2time(eph->week,data[27]),toc);

				eph->code=(int)data[20];      /* GPS: codes on L2 ch */
				eph->svh =(int)data[24];      /* sv health */
				eph->sva=uraindex(data[23]);  /* ura (m->index) */
				eph->flag=(int)data[22];      /* GPS: L2 P data flag */

				eph->tgd[0]=   data[25];      /* TGD */
				if (sys==SYS_GPS) {
					eph->fit=data[28];        /* fit interval (h) */
				}
				else {
					eph->fit=data[28]==0.0?1.0:2.0; /* fit interval (0:1h,1:>2h) */
				}
			}
			// Code for other satellite systems can be found in decode_eph in rinex.c in rtklib

			if (eph->iode<0||1023<eph->iode) {
				funcLogHelperM("rinex nav invalid: sat " << sat << " iode " << eph->iode);
			}
			if (eph->iodc<0||1023<eph->iodc) {
				funcLogHelperM("rinex nav invalid: sat " << sat << " iodc " << eph->iodc);
			}
			return 1;

		}



		/* satellite system+prn/slot number to satellite number ------------------------
		 * Taken from rtkcmn.c in rtklib
		 * convert satellite system+prn/slot number to satellite number
		 * args   : int    sys       I   satellite system (SYS_GPS,SYS_GLO,...)
		 *          int    prn       I   satellite prn/slot number
		 * return : satellite number (0:error)
		 *-----------------------------------------------------------------------------*/
		int satno(int sys, int prn)
		{
			if (prn<=0) return 0;
			switch (sys) {
				case SYS_GPS:
					if (prn<MINPRNGPS||MAXPRNGPS<prn) return 0;
					return prn-MINPRNGPS+1;
				/*
				case SYS_GLO:
					if (prn<MINPRNGLO||MAXPRNGLO<prn) return 0;
					return NSATGPS+prn-MINPRNGLO+1;
				case SYS_GAL:
					if (prn<MINPRNGAL||MAXPRNGAL<prn) return 0;
					return NSATGPS+NSATGLO+prn-MINPRNGAL+1;
				case SYS_QZS:
					if (prn<MINPRNQZS||MAXPRNQZS<prn) return 0;
					return NSATGPS+NSATGLO+NSATGAL+prn-MINPRNQZS+1;
				case SYS_CMP:
					if (prn<MINPRNCMP||MAXPRNCMP<prn) return 0;
					return NSATGPS+NSATGLO+NSATGAL+NSATQZS+prn-MINPRNCMP+1;
				case SYS_LEO:
					if (prn<MINPRNLEO||MAXPRNLEO<prn) return 0;
					return NSATGPS+NSATGLO+NSATGAL+NSATQZS+NSATCMP+prn-MINPRNLEO+1;
				case SYS_SBS:
					if (prn<MINPRNSBS||MAXPRNSBS<prn) return 0;
					return NSATGPS+NSATGLO+NSATGAL+NSATQZS+NSATCMP+NSATLEO+prn-MINPRNSBS+1;
					*/
				default:
					funcLogHelperM("satellite constellation " << sys << " unsupported");
					return -1;

			}
			return 0;
		}

		/* satellite number to satellite system ----------------------------------------
		 * Taken from rtkcmn.c in rtklib
		 * convert satellite number to satellite system
		 * args   : int    sat       I   satellite number (1-MAXSAT)
		 *          int    *prn      IO  satellite prn/slot number (NULL: no output)
		 * return : satellite system (SYS_GPS,SYS_GLO,...)
		 *-----------------------------------------------------------------------------*/
		int satsys(int sat, int *prn)
		{
			int sys=SYS_NONE;
			if (sat<=0||MAXSAT<sat) sat=0;
			else if (sat<=NSATGPS) {
				sys=SYS_GPS; sat+=MINPRNGPS-1;
			}
			/*
			else if ((sat-=NSATGPS)<=NSATGLO) {
				sys=SYS_GLO; sat+=MINPRNGLO-1;
			}
			else if ((sat-=NSATGLO)<=NSATGAL) {
				sys=SYS_GAL; sat+=MINPRNGAL-1;
			}
			else if ((sat-=NSATGAL)<=NSATQZS) {
				sys=SYS_QZS; sat+=MINPRNQZS-1;
			}
			else if ((sat-=NSATQZS)<=NSATCMP) {
				sys=SYS_CMP; sat+=MINPRNCMP-1;
			}
			else if ((sat-=NSATCMP)<=NSATLEO) {
				sys=SYS_LEO; sat+=MINPRNLEO-1;
			}
			else if ((sat-=NSATLEO)<=NSATSBS) {
				sys=SYS_SBS; sat+=MINPRNSBS-1;
			}*/
			else sat=0;
			if (prn) *prn=sat;
			return sys;
		}



		/* ura value (m) to ura index
		 *
		 * Taken from rinex.c
		 */
		static int uraindex(double value)
		{
		    int i;
		    for (i=0;i<15;i++) if (dsp::utils::URA_EPH[i]>=value) break;
		    return i;
		}

		/** set system mask
		 *
		 *  Taken from rinex.c
		 */
		int set_sysmask(const char *opt)
		{
		    const char *p;
		    int mask=SYS_NONE;

		    if (!(p=strstr(opt,"-SYS="))) return SYS_ALL;

		    for (p+=5;*p&&*p!=' ';p++) {
		        switch (*p) {
		            case 'G': mask|=SYS_GPS; break;
		            case 'R': mask|=SYS_GLO; break;
		            case 'E': mask|=SYS_GAL; break;
		            case 'J': mask|=SYS_QZS; break;
		            case 'C': mask|=SYS_CMP; break;
		            case 'I': mask|=SYS_IRN; break;
		            case 'S': mask|=SYS_SBS; break;
		        }
		    }
		    return mask;
		}



	}
}
}


