//optimazation parameteres for maf estimation
#define MAF_START 0.3
#define MAF_ITER 20

#include <string>
#include <iostream>

#include <cstdio>
// #include <cstdlib>
// #include <cstring>
#include <cmath>
// #include <limits>
#include <vector>
//#include <sys/stat.h>


double misTol=0.05;

double emFrequency(double *loglike,int numInds, int iter,double start,int *keep,int keepInd){

  if(keepInd == 0)
    return 0.0;
  
  float W0;
  float W1;
  float W2;
  // fprintf(stderr,"start=%f\n",start);
  float p=(float)start;
  float temp_p=(float)start;
  double accu=0.00001;
  double accu2=0;
  float sum;


  int it=0;
  
  for(it=0;it<iter;it++){
    sum=0;
    for(int i=0;i<numInds;i++){
      if(keep[i]==0)
        continue;
      W0=(loglike[i*3+0])*(1-p)*(1-p);
      W1=(loglike[i*3+1])*2*p*(1-p);
      W2=(loglike[i*3+2])*p*p;
      sum+=(W1+2*W2)/(2*(W0+W1+W2));
      //  fprintf(stderr,"%f %f %f\n",W0,W1,W2);
      if(0&&std::isnan(sum)){
	//fprintf(stderr,"PRE[%d]: W %f\t%f\t%f sum=%f\n",i,W0,W1,W2,sum);
	exit(0);
     }
    }

    p=sum/keepInd;
    // fprintf(stderr,"it=%d\tp=%f\tsum=%f\tkeepInd=%d\n",it,p,log(sum),keepInd);
    if((p-temp_p<accu&&temp_p-p<accu)||(p/temp_p<1+accu2&&p/temp_p>1-accu2))
      break;
    temp_p=p;
  }



  if(std::isnan(p)){
    fprintf(stderr,"[%s] caught nan will not exit\n",__FUNCTION__);
    fprintf(stderr,"logLike (3*nInd). nInd=%d\n",numInds);
    //print_array(stderr,loglike,3*numInds);
    fprintf(stderr,"keepList (nInd)\n");
    //print_array(stderr,keep,numInds);
    fprintf(stderr,"used logLike (3*length(keep))=%d\n",keepInd);

    for(int ii=0;0&&ii<numInds;ii++){
      if(keep!=NULL && keep[ii]==1)
	    fprintf(stderr,"1\t");
	for(int gg=0;gg<3;gg++)
	  fprintf(stderr,"%f\t",loglike[ii*3+gg]);
      fprintf(stderr,"\n");
    }
    sum=0;
    for(int i=0;i<numInds;i++){
      if(keep!=NULL && keep[i]==0)
        continue;
      W0=(loglike[i*3+0])*(1-p)*(1-p);
      W1=(loglike[i*3+1])*2*p*(1-p);
      W2=(loglike[i*3+2])*p*p;
      sum+=(W1+2*W2)/(2*(W0+W1+W2));
      //fprintf(stderr,"p=%f W %f\t%f\t%f sum=%f loglike: %f\n",p,W0,W1,W2,sum,exp(loglike[i*3+2])*pow(1-p,2));
    }
    p=-999;
    // exit(0);
  }
  
  return(p);
}
 

int main(int argc, char **argv) {
    double minMaf = atof(argv[1]);
    std::string delim = "\t";
    std::string row;
    std::string site_id;
    bool header=true;
    int nind=0;
    double *genos;
    int * keeps;
    int nsite = 0;
    int nexcluded = 0;
    while(std::getline(std::cin, row)){
        /* if (nsites %100000 == 0) */
        /*     std::cerr << "log:" << nsites << '\r' */
        if (header)
            std::cout << row << '\n';
        auto start = 0U; 
        auto end = row.find(delim);
        int counter = 0;
        int beagle_idx = 0;
        while (end != std::string::npos){

            if(header){
                nind++;        
            } else if ( counter == 0){
                site_id = row.substr(start, end - start);
            } else if (counter>2) {
                genos[beagle_idx] = atof(row.substr(start, end - start).c_str());
                beagle_idx++;
            }

            start = end + delim.length();
            end = row.find(delim, start);
                
            if(!header && genos[beagle_idx]<0){
                fprintf(stderr,"##ERROR: Likelihoods must be positive\n");
                fprintf(stderr, "##ERROR: see site %s\n", site_id.c_str());
                exit(0);
            }
            counter++;
        }
        if (header){
            nind++;
            nind = (nind-3) / 3;
            std::cerr << "##LOG:nind " << nind << std::endl;
            std::cerr << "##LOG:minmaf " << minMaf << std::endl;

            genos = new double[3*nind]; 
            keeps = new int[nind];    
            header=false;
            continue;
        } else {
            genos[beagle_idx] = atof(row.substr(start, end - start).c_str());
        }

        
        for(int i=0;i<nind;i++){
            double tmpS = 0.0;
            for(int g=0;g<3;g++)
                tmpS += genos[i*3+g];
            if(!(tmpS>0)){
                fprintf(stderr,"##ERROR: The sum of likelihoods for a genotypes must be positive\n");
                fprintf(stderr,"##ERROR: individual %d at site %d has sum %f\n",i, nsite,tmpS);
                for(int g=0;g<3;g++)
                    std::cerr << "##ERROR: " << g << " " << genos[i*3 + g] << '\n';
                exit(0);
            }
        }
        int nkeep=0;
        for(int i=0;i<nind;i++){
            double mmin=std::min(genos[i*3],std::min(genos[i*3+1],genos[i*3+2]));
            double mmax=std::max(genos[i*3],std::max(genos[i*3+1],genos[i*3+2]));
            if(fabs(mmax-mmin)<misTol){
                // skip site
                keeps[i] = 0;
                continue; 
            } else {
                keeps[i] = 1;
                nkeep++;
            }
        }
      
        double maf = emFrequency(genos,nind,MAF_ITER,MAF_START,keeps,nkeep);

        if (maf > minMaf && maf < 1-minMaf){
            std::cout << row << '\n';
        }else {
            fprintf(stderr, "SITEEXCLUDED %s maf: %f\n", site_id.c_str(), maf);
            nexcluded++;
            // excluded 
        }
    nsite++;
    }
    fprintf(stderr, "##LOG: parsed %d sites. Excluded %d sites with maf < %f\n", nsite, nexcluded, minMaf); 
    free(genos);
    free(keeps);
}
