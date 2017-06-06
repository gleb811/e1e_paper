#include <TMath.h>
#include <sstream>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <stdio.h> 
#include <math.h>
#include "mom_corr.h"
 
   Double_t mom_corr::correct_thel_e1_2039_2250_feb09(double pe,double theta,double phi)
    {
      //IMPLICIT none
      double phisec,fact,correct_thel_e1_2039_2250_feb09;
      double d,x,y,x2,y2,f1,f2;
      int nsector;
      phisec=(fmod(phi+30,60.)-30.);
     //phisec=((phi+30)%60.-30.);
      fact=1.0;
      nsector=0;
      if(phi>330 || phi <=30.) nsector=1;
      if(phi>30 && phi <= 90.) nsector=2;
      if(phi>90 && phi <= 150.) nsector=3;
      if(phi>150 && phi <= 210.) nsector=4;
      if(phi>210 && phi <= 270.) nsector=5;
      if(phi>270 && phi <= 330.) nsector=6;
//sector1
      if(nsector == 1){

        if ((theta >16.) && (theta < 32.5)){
             x=0.2051420E-01-.5001128E-03*(theta);
             y=-.2507647E-02+0.9725749E-04*(theta);
             x2=0.1526173E-03+0.4581319E-05*(theta)+0.1142428E-03*TMath::Sin(449./theta+8.054586);
             fact=x+y*(phisec)+x2*((phisec)*(phisec));
        }
        else if ((theta >32.5) && (theta < 44.)){
             x=6.252384-.4711988*(theta)+0.1178777E-01*((theta)*(theta))-.9834436E-04*((theta)*(theta)*(theta));
             y=-.2758633+0.2160061E-01*(theta)-.5579384E-03*((theta)*(theta))+0.4775051E-05*((theta)*(theta)*(theta));
             x2=-.3277342E-03+0.1174220E-04*(theta)-.2276303E-03*TMath::Sin(449./theta-8.186739);
             fact=x+y*(phisec)+x2*((phisec)*(phisec));
       }else fact=0;
     }
//sector2           
     if(nsector == 2){

        if ((theta >16.) && (theta < 32.5)){
            x=-.1339374+0.1013706E-01*(theta)-.1723397E-03*((theta)*(theta));
            y=-.1638717E-02+0.9508414E-04*(theta);
            x2=0.2754756E-03+0.7588073E-07*(theta)-.1080132E-03*TMath::Sin(449./theta-1.257899);
            fact=x+y*(phisec)+x2*((phisec)*(phisec));
        }
        else if ((theta>32.5) && (theta <44.)){
            x=5.225882-.4008103*(theta)+0.1023979E-01*((theta)*(theta))-.8721992E-04*((theta)*(theta)*(theta));
            y=-.6548676+0.5129098E-01*(theta)-.1329902E-02*((theta)*(theta))+0.1144865E-04*((theta)*(theta)*(theta));
            x2=-.3514004E-01+0.1304112E-02*(theta)+0.2458228E-04*((theta)*(theta))-.7083146E-06*((theta)*(theta)*(theta))-.2314476E-07*((theta)*(theta)*(theta)*(theta))+0.4713921E-09*((theta)*(theta)*(theta)*(theta)*(theta));
            fact=x+y*(phisec)+x2*((phisec)*(phisec));
       }else fact=0;
     }
//sector3
     if(nsector == 3){


        if ((theta >16.) && (theta <44.)){
           x=0.5706769-.9922977E-01*(theta)+0.5942143E-02*((theta)*(theta))-.1460487E-03*((theta)*(theta)*(theta))+0.1269023E-05*((theta)*(theta)*(theta)*(theta));
           y=0.7726650E-02-.1776088E-03*(theta);
           x2=0.5845903E-03-.1155695E-04*(theta)-.1278195E-03*TMath::Sin(449./theta-1.534029);
           fact=x+y*(phisec)+x2*((phisec)*(phisec));
        }else fact=0;



     }
//sector4

     if(nsector == 4){


        if ((theta>16.)&& (theta <44.)){
          x=1.767244-.2542930*(theta)+0.1246248E-01*((theta)*(theta))-.1430444E-03*((theta)*(theta)*(theta))-.6262982E-05*((theta)*(theta)*(theta)*(theta))+0.1960881E-06*((theta)*(theta)*(theta)*(theta)*(theta))-.1566946E-08*((theta)*(theta)*(theta)*(theta)*(theta)*(theta));
          y=-.1807822E-01+0.3758978E-02*(theta)-.2207941E-03*((theta)*(theta))+0.5200121E-05*((theta)*(theta)*(theta))-.4353816E-07*((theta)*(theta)*(theta)*(theta));
          x2=0.4989405E-03-.8992112E-05*(theta)+0.1226031E-03*TMath::Sin(449./theta+1.448175);
          fact=x+y*(phisec)+x2*((phisec)*(phisec));
        }else fact=0;


     }
//sector5
     if(nsector == 5){


        if ((theta >16.)&&(theta <32.5)){
          x=0.6389222E-01-.4801639E-03*(theta);
          y=0.3421679E-02-.8440115E-04*(theta);
          x2=0.3923497E-03-.4957983E-05*(theta)-.9631179E-04*TMath::Sin(449./theta-1.264236);
          fact=x+y*(phisec)+x2*((phisec)*(phisec));
        }
        else if ((theta >32.5) && (theta < 44.)){
          x=24.81553-1.917555*(theta)+0.4922545E-01*((theta)*(theta))-.4198686E-03*((theta)*(theta)*(theta));

          y=17.04936-1.342118*(theta)+0.2436030E-01*((theta)*(theta))+0.5353051E-03*((theta)*(theta)*(theta))-.2185274E-04*((theta)*(theta)*(theta)*(theta))+0.1889016E-06*((theta)*(theta)*(theta)*(theta)*(theta));
          x2=0.2119622E-02-.5375929E-04*(theta)-.2181619E-03*TMath::Sin(449./theta-.1648105);
          fact=x+y*(phisec)+x2*((phisec)*(phisec));

        }else fact=0;



    }
//sector6
    if(nsector == 6){

       if ((theta >16.) && (theta <32.5)){
          x=0.2883536E-01-.4278146E-04*(theta);
          y=0.4266386E-02-.7680542E-04*(theta);
          x2=0.2714948E-03+0.6075697E-06*(theta)+0.7332848E-04*TMath::Sin(449./theta+8.132975);
          fact=x+y*(phisec)+x2*((phisec)*(phisec));
        }
        else if ((theta > 32.5)&& (theta <44.)){
          x=8.539058-.6286885*(theta)+0.1536358E-01*((theta)*(theta))-.1245493E-03*((theta)*(theta)*(theta));

          y=-.1018162E-02+0.6171986E-04*(theta)+0.8155919E-06*((theta)*(theta))+0.1453145E-07*((theta)*(theta)*(theta))+0.2333810E-09*((theta)*(theta)*(theta)*(theta))-.2149097E-10*((theta)*(theta)*(theta)*(theta)*(theta));
          x2=0.1190329E-02-.2835207E-04*(theta)-.1359806E-03*TMath::Sin(449./theta-.8999598);
          fact=x+y*(phisec)+x2*((phisec)*(phisec));
        }else fact=0;
    }

   correct_thel_e1_2039_2250_feb09=theta-fact;
   return correct_thel_e1_2039_2250_feb09;

  }
Double_t mom_corr::correct_pel_e1_2039_2250_feb09(double pe,double theta,double phi){
      //IMPLICIT none
      double phisec,fact,correct_pel_e1_2039_2250_feb09;
      double d,x,y,x2,y2,f1,f2;
      int nsector;
      phisec=(fmod(phi+30,60.)-30.);
      fact=1.0;
      nsector=0;
      if(phi>330 || phi <=30.) nsector=1;
      if(phi>30 && phi <= 90.) nsector=2;
      if(phi>90 && phi <= 150.) nsector=3;
      if(phi>150 && phi <= 210.) nsector=4;
      if(phi>210 && phi <= 270.) nsector=5;
      if(phi>270 && phi <= 330.) nsector=6;

      if(nsector == 1){

       if ((theta >16.)&& (theta < 52.)){


          x=1.008854-.3422894E-03*(theta)-.6954598E-06*((theta)*(theta))+0.8156167E-07*((theta)*(theta)*(theta));
          y=-.1025022E-02+0.8506941E-05*(theta)+0.1281147E-05*((theta)*(theta))-.2162752E-07*((theta)*(theta)*(theta));

          x2=-.9028474E-04+0.8689741E-05*(theta)-.2580860E-06*((theta)*(theta))+0.2403976E-08*((theta)*(theta)*(theta));
          fact=(x+y*(phisec))+x2*((phisec)*(phisec));
       }else fact=1.;
      }
      if(nsector == 2){

       if ((theta >16.)&& (theta < 52.)){

          x=0.9984451+0.5805102E-04*(theta)-.2909986E-05*((theta)*(theta))+0.4793942E-07*((theta)*(theta)*(theta));

          y=-.1578927E-02+0.9188605E-04*(theta)-.1745132E-05*((theta)*(theta))+0.8876872E-08*((theta)*(theta)*(theta));
          x2=-.4820524E-04+0.5224319E-05*(theta)-.1622990E-06*((theta)*(theta))+0.1537213E-08*((theta)*(theta)*(theta));
          fact=(x+y*(phisec))+x2*((phisec)*(phisec));
       }else fact=1.;
      }
     if(nsector == 3){

       if ((theta >16.)&& (theta < 52.)){


          x=1.017301-.2568450E-02*(theta)+0.1060422E-03*((theta)*(theta))-.1080502E-05*((theta)*(theta)*(theta))-.1366504E-07*((theta)*(theta)*(theta)*(theta))+0.2249102E-09*((theta)*(theta)*(theta)*(theta)*(theta));

          y=-.2755583E-02+0.2811991E-03*(theta)-.6515745E-05*((theta)*(theta))-.5426127E-07*((theta)*(theta)*(theta))+0.3184020E-08*((theta)*(theta)*(theta)*(theta))-.2669610E-10*((theta)*(theta)*(theta)*(theta)*(theta));

          x2=-.8553622E-04+0.8231474E-05*(theta)-.2452097E-06*((theta)*(theta))+0.2279487E-08*((theta)*(theta)*(theta));
          fact=(x+y*(phisec))+x2*((phisec)*(phisec));
       }else fact=1.;
    }
     if(nsector == 4){

       if ((theta >16.)&& (theta < 52.)){

          x=1.016233-.7399115E-03*(theta)+0.8491787E-05*((theta)*(theta));

          y=-.3036835E-02+0.3991334E-03*(theta)-.1802434E-04*((theta)*(theta))+0.3370880E-06*((theta)*(theta)*(theta))-.2288362E-08*((theta)*(theta)*(theta)*(theta));
          x2=-.4231546E-04+0.4222708E-05*(theta)-.1275686E-06*((theta)*(theta))+0.1193804E-08*((theta)*(theta)*(theta));
          fact=(x+y*(phisec))+x2*((phisec)*(phisec));
       }else fact=1.;
    }
     if(nsector == 5){

       if ((theta >16.)&& (theta < 52.)){


          x=0.9879246+0.1914278E-02*(theta)-.7166770E-04*((theta)*(theta))+0.8048290E-06*((theta)*(theta)*(theta));
          y=-.1823859E-02+0.1565340E-03*(theta)-.4768085E-05*((theta)*(theta))+0.4759397E-07*((theta)*(theta)*(theta));

          x2=-.5335819E-04+0.5686952E-05*(theta)-.1885548E-06*((theta)*(theta))+0.1940291E-08*((theta)*(theta)*(theta));
          y2= -.2186460E-07-.6860175E-09*(theta);
          d=0.5041720E-08+0.1139555E-09*(theta);
          fact=(x+y*(phisec))+x2*((phisec)*(phisec))+y2*((phisec)*(phisec)*(phisec))+d*((phisec)*(phisec)*(phisec)*(phisec));
       }else fact=1.;
     }
    if(nsector == 6){

       if ((theta >16.)&& (theta < 52.)){

          x=1.007478-.2812020E-03*(theta)+0.3874791E-05*((theta)*(theta));

          y=-.7835627E-03+0.9022003E-04*(theta)-.2670348E-05*((theta)*(theta))+0.2166390E-07*((theta)*(theta)*(theta));
          x2=0.3509918E-05-.1477360E-06*(theta);
          fact=(x+y*(phisec))+x2*((phisec)*(phisec));
       }else fact=1.;
    }

      correct_pel_e1_2039_2250_feb09=fact*pe;
      return correct_pel_e1_2039_2250_feb09;
   }
