utils::globalVariables(c("cfm1", "cfm2","cfm3","cfp1","cfp2","cfp3","a0","a1","a2","a3","wls","dbloss","size","propa","a4"))
#' @import stats

genProps<-function(){
#shared standard params
  prop<-list()
  prop$reference_att = NA;
  prop$dist = NA;
  prop$antenna_height = 0;
  prop$wave_number = NA;
  prop$terrain_irregularity = NA;
  prop$surface_refractivity = NA;
  prop$earth_eff_curvature = NA;
  prop$surface_xfr_impedance = NA;
  prop$effective_height= 0;
  prop$horizon_dist = 0;
  prop$horizon_elev_angle = 0;
  prop$error_id = NA;
  prop$errmsg = NA;
  prop$mode = NA;

  #shared scattering params
  propa<-list()
  propa$los_distance = NA;
  propa$scatter_dist = NA;

  #line of sight coefficients
  propa$ael = NA;
  propa$ak1 = NA;
  propa$ak2 = NA;

  #diffraction coefficients
  propa$aed = NA;
  propa$emd = NA;

  #scatter coefficients
  propa$aes = NA;
  propa$ems = NA;

  propa$smooth_horizon_dist = 0;
  propa$total_horizon_dist = NA;
  propa$total_bending_angle = NA;

  #shared variability parameters
  propv<-list()
  propv$sgc_confidence = NA;
  propv$ctrl_state = NA;
  propv$mode_variability = NA;
  propv$climate_id = NA;
  returnme<-list(prop,propa,propv)
  names(returnme)<-c("prop","propa","propv")
  returnme
}

fKnifeAtt<-function(v2){ #function output =
  #v2 is v^2, approximate Fresnel integral in dB using the following formula
  if(v2<5.76){
    output=6.02+9.11*sqrt(v2)-1.27*v2;
  }else{
    output=12.953+4.343*log(v2);
  }
  output
}

fHeightGain<-function( x,  pk){ #function output =
  #This is an approximation for diffraction around a smooth earth using x as height
  if(x<200.0){
    w=-log(pk);
    if (( pk < 1e-5) | ((x*w^3) > 5495.0 )){
      output=-117.0;
      if (x>1.0){
        output=17.372*log(x)+output;
      }
    }else{
      output=2.5e-5*x^2/pk+20*log(pk)-15;
    }
  }else{
    output=0.05751*x-4.343*log(x);
    if(x<2000.0){
      w=0.0134*x*exp(-0.005*x);
      output=(1.0-w)*output+w*(17.372*log(x)-117.0);
    }

  }
  output
}

fDiffractionAtt<-function(d, prop, propa, tempvars=list()){ #function [output, prop, propa] =
  #Diffraction Attenuation at distance D. Uses a convex combination
  #of smooth earth diffraction and double knife-edge diffraction. an
  #initial call with d=0 sets up initial constants

  #Prepare initial diffraction constants (Section 11)
  if(d==0){
    q=prop$antenna_height[1]*prop$antenna_height[2];
    qk=prop$effective_height[1]*prop$effective_height[2]-q;
    if(prop$mode<0){#point to point
      q=q+10;
    }

    wd1 = sqrt(1+qk/q);
    xd1 = propa$total_horizon_dist+propa$total_bending_angle/prop$earth_eff_curvature;
    q = (1-0.8*exp(-propa$los_distance/50e3))*prop$terrain_irregularity;
    q = q*0.78*exp(-(q/16)^0.25);
    afo = min(15,2.171*log(1+4.77e-4*prop$antenna_height[1]*prop$antenna_height[2] * prop$wave_number*q));
    qk = 1/abs(prop$surface_xfr_impedance);
    aht = 20;
    xht = 0;
    for (j in 1:2){
      a=0.5*prop$horizon_dist[j]^2/prop$effective_height[j];
      wa=(a*prop$wave_number)^(1/3);
      pk=qk/wa;  #modify for static var
      q=(1.607-pk)*151*wa*prop$horizon_dist[j]/a;
      xht =xht + q;
      aht =aht + fHeightGain(q,pk);
      #output=0.0;

      #set output
      output=list()
      output$wd1 = wd1;
      output$xd1 = xd1;
      output$afo = afo;
      output$qk = qk;
      output$aht = aht;
      output$xht = xht;
    }

    returnme<-list(output, prop, propa)
    names(returnme)<-c("tempvars", "prop", "propa")
    return(returnme)

  }else{
    #Calculate diffraction Attenuation
    th = propa$total_bending_angle+d*prop$earth_eff_curvature;
    ds = d-propa$total_horizon_dist;
    q = 0.0795775*prop$wave_number*ds*th^2;
    output = fKnifeAtt(q*prop$horizon_dist[1]/(ds+prop$horizon_dist[1]))+fKnifeAtt(q*prop$horizon_dist[2]/(ds+prop$horizon_dist[2]));
    a = ds/th;
    wa = (a*prop$wave_number)^(1/3);
    pk = tempvars$qk/wa;
    q = (1.607-pk)*151.0*wa*th + tempvars$xht;
    ar = 0.05751*q-4.343*log(q) - tempvars$aht;
    q = (tempvars$wd1+tempvars$xd1/d) * min(((1.0-0.8*exp(-d/50e3))*prop$terrain_irregularity*prop$wave_number), 6283.2 );
    wd = 25.1/(25.1+sqrt(q));
    output = ar*wd+(1-wd)*output+tempvars$afo;

    returnme<-list(output, prop, propa,tempvars)
    names(returnme)<-c("output", "prop", "propa","tempvars")
    return(returnme)

  }

}

fLongleyRicePropSetup<-function( fmhz, zsys, en0, ipol, eps, sgm, prop){
  gma = 157e-9; #fixed for Earth.
  prop$wave_number = fmhz/47.7; #assums speed of light propagation speed
  prop$surface_refractivity = en0;
  if(zsys!=0.0) prop$surface_refractivity = prop$surface_refractivity*exp(-zsys/9460);

  prop$earth_eff_curvature = gma*(1-0.04665*exp(prop$surface_refractivity/179.3));

  zq = eps + 1i*376.73*sgm/prop$wave_number; #changed from 376.62 to be more accurate and use Constants_.m
  prop$surface_xfr_impedance = sqrt(zq-1.0);
  if(ipol!=0.0) prop$surface_xfr_impedance = prop$surface_xfr_impedance/zq;
  prop
}

fCalcTerrainIrregularity<-function(path_profile, x1, x2){ #d1thxv

  np=(path_profile[1]);           #was int, will be int, since its the number of points
  xa=x1/path_profile[2];          #
  xb=x2/path_profile[2];
  d1thxv=0;
  if((xb-xa)<2.0){  # exit out
    return(d1thxv);
  }
  ka=floor(0.1*(xb-xa+8));       #was int
  ka=min(max(4,ka),25);
  n=10*ka-5;
  kb=n-ka+1;
  sn=n-1;

  #s = zeros(1,n+2);
  s <- rep(0,n+2)
  s[1]=sn;
  s[2]=1;
  xb=(xb-xa)/sn;
  k=floor(xa+1.0);                 #was int
  xa=xa-k;
  for (j in 1:(n)) {                    #MOD FROM j=0 to n-1

    while(xa>0 & k<np){
      xa=xa-1;
      k = k+1;
    }

    s[j+2]=path_profile[k+3]+(path_profile[k+3]-path_profile[k+2])*xa;
    xa=xa+xb;

  }

  templist<-fLsqFit(s,0.0,sn); #[xa,xb]=
  names(templist)<-c("xa","xb")
  list2env(templist,envir=environment())
  xb=(xb-xa)/sn;
  for (j in 1:(n)){         #mod from J=0 to n-1
    s[j+2]=s[j+2]-xa;
    xa=xa+xb;
}

  d1thxv = quantile(s(3:size(s,2)),1-(ka/n)) - quantile(s(3:size(s,2)),1-(kb/n)); #use matlab's quantile function
  d1thxv = d1thxv/(1-0.8*exp(-(x2-x1)/50e3));
  return(d1thxv)
}

fLongleyRiceArea<-function(kst, klimx, mdvarx,prop, propv){

  for (j in 1:2){
    if(kst[j]<=0){
      prop$effective_height[j]=prop$antenna_height[j];
    }else{
      q=4.0;
      if (kst[j]!=1){
        q=9.0;
      }
      if (prop$antenna_height[j]<5.0){
        q=q*sin(pi/10*prop$antenna_height[j]);
      }
      prop$effective_height[j]=prop$antenna_height[j]+(1.0+q)*exp(-min(20.0,2.0*prop$antenna_height[j]/max(1e-3,prop$terrain_irregularity)));
    }
    q = sqrt(2.0*prop$effective_height[j]/prop$earth_eff_curvature);
    prop$horizon_dist[j] = q * exp(-0.07*sqrt(prop$terrain_irregularity/max(prop$effective_height[j],5)));
    prop$horizon_elev_angle[j] = (0.65*prop$terrain_irregularity * (q/prop$horizon_dist[j]-1) - 2*prop$effective_height[j]) /q;
  }


  prop$mode=1; #area mode
  propv$ctrl_state=max(propv$ctrl_state,3);
  if(mdvarx>=0){
    propv$mode_variability=mdvarx;
    propv$ctrl_state=max(propv$ctrl_state,4);
  }
  if(klimx>0){
    propv$climate_id=klimx;
    propv$ctrl_state=5;
  }
  returnme<-list(prop,propv)
  names(returnme)<-c("prop","propv")
  returnme
}

fLosAtt<-function(d, prop, propa, wls ){ #function [alosv, prop, propa] =
  if(d==0.0){
    #Generate initial constants
    wls=0.021/(0.021+prop$wave_number*prop$terrain_irregularity/max(10e3,propa$los_distance));
    #alosv=0.0;
    alosv = wls;
  }else{
    q=(1-0.8*exp(-d/50e3))*prop$terrain_irregularity;
    s=0.78*q*exp(-(q/16.0)^0.25);
    q=prop$effective_height[1]+prop$effective_height[2];
    sps=q/sqrt(d^2+q^2);
    r=(sps-prop$surface_xfr_impedance)/(sps+prop$surface_xfr_impedance)*exp(-min(10,prop$wave_number*s*sps));
    q=abs(r)^2;
    if (q<0.25 | q<sps){
      r=r*sqrt(sps/q);
    }
    alosv=propa$emd*d+propa$aed;
    q=prop$wave_number*prop$effective_height[1]*prop$effective_height[2]*2/d;
    if  (q > pi/2){
      q=pi-2.4649/q;
    }
    #alosv=(-4.343*log(obj.abq_alos((cos(q)-1i*sin(q))+r))-alosv) * obj.wls+alosv;
    alosv =(-4.343*log(abs(cos(q)-1i*sin(q)+r)^2)-alosv)*wls + alosv;
  }
  returnme<-list(alosv, prop, propa)
  names(returnme)<-c("alosv", "prop", "propa")
  return(returnme)
}

fH01<-function(r, et){ #function h0fv =
  a=c(25, 80, 177, 395, 705);
  b=c(24, 45,  68,  80, 105);
  it = floor(et);
  if(it<=0){
    it=1;
    q=0.0;

  }else{
    if(it>=5){
      it=5;
      q=0.0;

    }else{
      q=et-it;
    }
  }

  x= r^(-2);

  h0fv=4.343*log((a[it]*x+b[it])*x+1.0);
  if(q!=0.0){
    h0fv=(1.0-q)*h0fv+q*4.343*log((a[it+1]*x+b[it+1])*x+1.0);
  }
  h0fv
}

fFtheta_d<-function(td){ #function output
  a = c(   133.4,    104.6,     71.8);
  b = c(0.332e-3, 0.212e-3, 0.157e-3);
  c = c(  -4.343,   -1.086,    2.171);
  if(td<=10e3){
    i=1;
  }else if(td<=70e3){
    i=2;
  }else{
    i=3;
  }

    output =  a[i]+b[i]*td+c[i]*log(td);
    output
  }

fTroposcatterAtt<-function(d, prop, propa, tempvars){ #function [ascatv, prop, propa] =
  ascatv<-list()
  if(d==0){ # set up initial coeficents
    ascatv$ad=prop$horizon_dist[1]-prop$horizon_dist[2];
    ascatv$rr=prop$effective_height[2]/prop$effective_height[1];
    if (ascatv$ad<0){
      ascatv$ad = -ascatv$ad;
      ascatv$rr = 1/ascatv$rr;
    }
    ascatv$etq = (5.67e-6*prop$surface_refractivity-2.32e-3)*prop$surface_refractivity+0.031;
    ascatv$h0s = -15;
    #ascatv = 0;
    return(ascatv)
  }else{
    if(tempvars$h0s>15.0){
      h0 = tempvars$h0s;

    }else{
      th=prop$horizon_elev_angle[1]+prop$horizon_elev_angle[2]+d*prop$earth_eff_curvature;
      r2=2*prop$wave_number*th;
      r1=r2*prop$effective_height[1];
      r2=r2*prop$effective_height[2];
      if ((r1<0.2) & (r2<0.2)){
        ascatv =  1001;  # <==== early return - function is undefined
        return(ascatv);
      }

      ss = (d-tempvars$ad)/(d+tempvars$ad);
      q = tempvars$rr/ss;
      ss = max(0.1,ss);
      q = min(max(0.1,q),10);
      z0 = (d-tempvars$ad)*(d+tempvars$ad)*th*0.25/d;
      et = (tempvars$etq*exp(-min(1.7,z0/8e3)^6)+1)*z0/1.7556e3;
      ett = max(et,1);
      h0 = (fH01(r1,ett)+fH01(r2,ett))*0.5;  #function H01 as used in alg 6.13
      h0 = h0 + min(h0,(1.38-log(ett))*log(ss)*log(q)*0.49);
      h0 = max(h0,0);
      if(et<1){
        h0=et*h0+(1-et)*4.343*log(((1+sqrt(2)/r1) * (1+sqrt(2)/r2))^2*(r1+r2)/(r1+r2+2*sqrt(2)));
      }

      if(h0>15 & tempvars$h0s>=0){
        h0=tempvars$h0s;
      }
    }
    tempvars$h0s=h0;
    th=propa$total_bending_angle+d*prop$earth_eff_curvature;
    ascatv=fFtheta_d(th*d)+4.343*log(47.7*prop$wave_number*(th)^4.0)-0.1 * (prop$surface_refractivity-301.0)*exp(-th*d/40e3)+h0; #fFtheta_d is scattered fields
    returnme<-list(ascatv,prop,propa)
    names(returnme)<-c("ascatv","prop","propa")
    returnme
  }

}

fLongleyRiceProp<-function(d, prop, propa){ #[prop, propa]
  #Main Longley Rice Propagation program, implements Longley-Rice model.
  LR_vars<-list() #needed?
  if(prop$mode!=0){ #mode is program flow
    #point to point, or area beginning.

    propa$smooth_horizon_dist[1:2]=sqrt(2*prop$effective_height[c(1:2)]/prop$earth_eff_curvature);

    propa$los_distance=propa$smooth_horizon_dist[1]+propa$smooth_horizon_dist[2];
    propa$total_horizon_dist=prop$horizon_dist[1]+prop$horizon_dist[2];
    propa$total_bending_angle=max(prop$horizon_elev_angle[1]+prop$horizon_elev_angle[2],-propa$total_horizon_dist*prop$earth_eff_curvature);
    LR_vars$wlos = 0;
    LR_vars$wscat = 0;

    #check parameters
    if(prop$wave_number<0.838){
      prop$error_id=max(prop$error_id,1);
      msg = paste0('Warning: Input Frequency ', (prop$wave_number*47.7), ' MHz is nearly out of range, for best results use frequency between 40 MHz and 10 GHz');
      prop$errmsg = rbind(prop$errmsg,msg);
    }
    if (prop$wave_number>210){
      prop$error_id=max(prop$error_id,1);
      msg = paste0('Warning: Input Frequency ', (prop$wave_number*47.7), ' MHz is nearly out of range, for best results use frequency between 40 MHz and 10 GHz');
      prop$errmsg = rbind(prop$errmsg,msg);
    }

    for (j in 1:2){
      if (prop$antenna_height[j]<1){
        prop$error_id=max(prop$error_id,1);
        msg = ('Warning: Antenna height less than 1m - Antennas should be between 1 and 1000m above ground\n');
        prop$errmsg = rbind(prop$errmsg,msg);
      }
      if (prop$antenna_height[j]>1000){
        prop$error_id=max(prop$error_id,1);
        msg = ('Warning: Antenna above 1000m, Antennas should be between 1 and 1000m above ground\n');
        prop$errmsg = rbind(prop$errmsg,msg);
      }
    }
    for (j in 1:2){
      if( abs(prop$horizon_elev_angle[j]) >200e-3){
        prop$error_id=max(prop$error_id,3);
        msg = ('Warning: Horizon elevation angle above 11.5 degrees, results probably invalid\n');
        prop$errmsg = rbind(prop$errmsg,msg);
      }
      if (prop$horizon_dist[j]<0.1*propa$smooth_horizon_dist[j]){
        prop$error_id=max(prop$error_id,3);
        msg = ('Warning: Horizon distance too short - Probably too close to a large hill. \n');
        prop$errmsg = rbind(prop$errmsg,msg);
      }
      if (prop$horizon_dist[j]>3*propa$smooth_horizon_dist[j]){
        prop$error_id=max(prop$error_id,3);
        msg = ('Warning: Horizon distance too long - Horizon is farther than smooth Earth horizon\n');
        prop$errmsg = rbind(prop$errmsg,msg);
      }
    }

    if (prop$surface_refractivity < 250){
      prop$error_id=4;
      msg = ('Warning: Surface refractivity less than 250, out of range');
      prop$errmsg = rbind(prop$errmsg,msg);
    }
    if (prop$surface_refractivity > 400){
      prop$error_id=4;
      msg = ('Warning: Surface refractivity over 400, out of range');
      prop$errmsg = rbind(prop$errmsg,msg);
    }
    if (prop$earth_eff_curvature < 75e-9){
      prop$error_id=4;
      msg = ("Warning: Earth's curvature less than 75e-9, out of range");
      prop$errmsg = rbind(prop$errmsg,msg);
    }
    if (prop$earth_eff_curvature > 250e-9){
      prop$error_id=4;
      msg = ("Warning: Earth's curvature over 250e-9, out of range");
      prop$errmsg = rbind(prop$errmsg,msg);
    }
    if (Re(prop$surface_xfr_impedance) <= abs(Im(prop$surface_xfr_impedance))){
      prop$error_id=4;
      msg = ('Warning: Surface transfer impedance real value less than imaginary value, out of range');
      prop$errmsg = rbind(prop$errmsg,msg);
    }
    if (prop$wave_number  < 0.419){
      prop$error_id=4;
      msg = ('Warning: Frequency less than 20Mhz, Frequency must be between 20MHz and 20 GHz');
      prop$errmsg = rbind(prop$errmsg,msg);
    }
    if (prop$wave_number  > 420){
      prop$error_id=4;
      msg = ('Warning: Frequency over 20GHz, Frequency must be between 20MHz and 20 GHz');
      prop$errmsg = rbind(prop$errmsg,msg);
    }


    for (j in 1:2){
      if(prop$antenna_height[j]<0.5){
        msg = ('Warning: Antenna height too low - Antennas must be between 0.5 and 3000m above ground\n');
        prop$errmsg = rbind(prop$errmsg,msg);
        prop$error_id=4;
      }
      if (prop$antenna_height[j]>3000){
        msg = ('Warning: Antenna height too high - Antennas must be between 0.5 and 3000m above ground\n');
        prop$errmsg = rbind(prop$errmsg,msg);
        prop$error_id=4;
      }

    }

    LR_vars$dmin=abs(prop$effective_height[1]-prop$effective_height[2])/200e-3;

    #Diffraction Coefficients (Section 9)
    list2env(fDiffractionAtt(0,prop,propa),envir = environment()) #[tempvars, prop, propa]
    LR_vars$xae=(prop$wave_number*prop$earth_eff_curvature^2)^-(1/3);
    d3=max(propa$los_distance,1.3787*LR_vars$xae+propa$total_horizon_dist);
    d4=d3+2.7574*LR_vars$xae;
    templist<-(fDiffractionAtt(d3, prop, propa, tempvars)) #[a3, prop, propa]
    names(templist)[1]<-"a3"
    list2env(templist,envir = environment())
    templist<-fDiffractionAtt(d4, prop, propa, tempvars) #[a4, prop, propa]
    names(templist)[1]<-"a4"
    list2env(templist,envir = environment())

    propa$emd=(a4-a3)/(d4-d3);

    propa$aed=a3-propa$emd*d3;

  }
  if (prop$mode>=0){
    #area mode, set flag for next round.
    prop$mode=0;
    prop$dist=d;
  }
  if(prop$dist>0){
    #ensure distance is positive
    if(prop$dist>1000e3){
      prop$error_id=max(prop$error_id,1);
      prop$errmsg = rbind(prop$errmsg, 'Warning: Distance over 1000km');
    }
    if(prop$dist<LR_vars$dmin){
      prop$error_id=max(prop$error_id,3);
      prop$errmsg = rbind(prop$errmsg, paste0('Warning: Distance less than #1.2fm, Distance is out of range -- ',LR_vars$dmin));
    }
    if(prop$dist<1e3 | prop$dist>2000e3){
      prop$error_id=4;
      prop$errmsg = rbind(prop$errmsg, 'Warning: Distance is out of range. Acceptable parameters are between 1km and 2000km')
    }
  }
  if(prop$dist<propa$los_distance){
    #We have line of sight, calc based on LOS
    if(!LR_vars$wlos){#This is to avoid recalculating if in area mode

      templist<-fLosAtt(0, prop, propa) #[wls, prop, propa]
      names(templist)[1]<-"wls"
      list2env(templist,envir = environment())
      d2=propa$los_distance;
      a2=propa$aed+d2*propa$emd;
      d0=1.908*prop$wave_number*prop$effective_height[1]*prop$effective_height[2];
      if(propa$aed>=0.0){
        d0=min(d0,0.5*propa$total_horizon_dist);
        d1=d0+0.25*(propa$total_horizon_dist-d0);

      }else{ d1=max(-propa$aed/propa$emd,0.25*propa$total_horizon_dist);
      }


      templist<-fLosAtt(d1, prop, propa, wls); #[a1, prop, propa]
      names(templist)[1]<-"a1"
      list2env(templist,envir = environment())
      templist<-NULL

      wq=0;
      if (d0<d1){

        templist<-fLosAtt(d0, prop, propa, wls) #[a0, prop, propa]
        names(templist)[1]<-"a0"
        list2env(templist,envir = environment())
        templist<-NULL

        q=log(d2/d0);
        propa$ak2=max(0.0,((d2-d0)*(a1-a0)-(d1-d0)*(a2-a0)) /((d2-d0)*log(d1/d0)-(d1-d0)*q));

        wq=propa$aed>=0.0 | propa$ak2>0.0;


        if(wq){

          propa$ak1=(a2-a0-propa$ak2*q)/(d2-d0);
          if(propa$ak1<0.0){
            propa$ak1=0.0;
            propa$ak2=max(a2-a0,0)/q;
            if(propa$ak2==0.0){
              propa$ak1=propa$emd;
            }
          }
        }
      }

      if(!wq){
        propa$ak1=max(a2-a1,0)/(d2-d1);
        propa$ak2=0.0;
        if(propa$ak1==0.0){
          propa$ak1=propa$emd;
        }
      }
      propa$ael=a2-propa$ak1*d2-propa$ak2*log(d2);
      LR_vars$wlos=1;

    }
    if(prop$dist>0.0){
      prop$reference_att=propa$ael+propa$ak1*prop$dist + propa$ak2*log(prop$dist);
    }
  }
  if(prop$dist<=0.0 | prop$dist>=propa$los_distance){
    #Troposcatter dominiant - calc based scattering
    if (!LR_vars$wscat){

      tempvars =fTroposcatterAtt(0,prop,propa);
      d5 = propa$total_horizon_dist+200e3;
      d6 = d5+200e3;
      a6 = fTroposcatterAtt(d6,prop,propa,tempvars)[[1]];
      a5 = fTroposcatterAtt(d5,prop,propa,tempvars)[[1]];
      if(a5<1000){
        propa$ems=(a6-a5)/200e3;
        propa$scatter_dist=max(propa$los_distance,max(propa$total_horizon_dist+0.3*LR_vars$xae * log(47.7*prop$wave_number),(a5-propa$aed-propa$ems*d5) /(propa$emd-propa$ems)));
        propa$aes=(propa$emd-propa$ems)*propa$scatter_dist+propa$aed;

      }else{
        propa$ems=propa$emd;
        propa$aes=propa$aed;
        propa$scatter_dist=10e6;
      }
      LR_vars$wscat=1;
    }
    if(prop$dist>propa$scatter_dist){
      prop$reference_att=propa$aes+propa$ems*prop$dist;
    }else{
      prop$reference_att=propa$aed+propa$emd*prop$dist;
    }

  }
  prop$reference_att=max(prop$reference_att,0.0);
  returnme<-list(prop,propa)
  names(returnme)<-c("prop","propa")
  returnme
}

fCurve<-function(c1, c2, x1, x2, x3, de){
    (c1+c2/(1.0+(((de-x2)/x3)^2.0)))*((de/x1)^2) /(1.0+((de/x1)^2))
}

fAttVariability<-function(zzt, zzl, zzc, prop, propv){ #function [avarv, prop, propv] =
  #computes the attenuation vs statistics (Section 28
  #When in area predicition mode one needs threefold quantile of
  #attenuation which corresponds to the fraction qt of time, ql of
  #locations, and qs as situations.
  #in point to point mode one needs only qt and qs, in this implementation for pt 2 pt, ql is set to .5.
  #Written to reduce duplicate calculations. on first entering set
  #ctrl_state  = 5. this will initialize all parameters
  #ctrl_state  = 4. use when modvar is changed
  #ctrl_state  = 3. use when freq is changed
  #ctrl_state  = 2. use when antenna heights are changed
  #ctrl_state  = 1. use when distances are changed.
  #the higher the lvar value, the more calculations to be performed.

  #climatic constants  (Section 29)
  #   equitorial, continental subtropical, maritime subtropical,
  #   desert, continental temperate, maritime over land, maritime
  #   over sea.
  bv1=c(-9.67,-0.62,1.26,-9.21,-0.62,-0.39,3.15);
  bv2=c(12.7,9.19,15.5,9.05,9.19,2.86,857.9);
  xv1=c(144.9e3,228.9e3,262.6e3,84.1e3,228.9e3,141.7e3,2222.e3);
  xv2=c(190.3e3,205.2e3,185.2e3,101.1e3,205.2e3,315.9e3,164.8e3);
  xv3=c(133.8e3,143.6e3,99.8e3,98.6e3,143.6e3,167.4e3,116.3e3);
  bsm1=c(2.13,2.66,6.11,1.98,2.68,6.86,8.51);
  bsm2=c(159.5,7.67,6.65,13.11,7.16,10.38,169.8);
  xsm1=c(762.2e3,100.4e3,138.2e3,139.1e3,93.7e3,187.8e3,609.8e3);
  xsm2=c(123.6e3,172.5e3,242.2e3,132.7e3,186.8e3,169.6e3,119.9e3);
  xsm3=c(94.5e3,136.4e3,178.6e3,193.5e3,133.5e3,108.9e3,106.6e3);
  bsp1=c(2.11,6.87,10.08,3.68,4.75,8.58,8.43);
  bsp2=c(102.3,15.53,9.60,159.3,8.12,13.97,8.19);
  xsp1=c(636.9e3,138.7e3,165.3e3,464.4e3,93.2e3,216.0e3,136.2e3);
  xsp2=c(134.8e3,143.7e3,225.7e3,93.1e3,135.9e3,152.0e3,188.5e3);
  xsp3=c(95.6e3,98.6e3,129.7e3,94.2e3,113.4e3,122.7e3,122.9e3);
  bsd1=c(1.224,0.801,1.380,1.000,1.224,1.518,1.518);
  bzd1=c(1.282,2.161,1.282,20.,1.282,1.282,1.282);
  bfm1=c(1.0,1.0,1.0,1.0,0.92,1.0,1.0);
  bfm2=c(0.0,0.0,0.0,0.0,0.25,0.0,0.0);
  bfm3=c(0.0,0.0,0.0,0.0,1.77,0.0,0.0);
  bfp1=c(1.0,0.93,1.0,0.93,0.93,1.0,1.0);
  bfp2=c(0.0,0.31,0.0,0.19,0.31,0.0,0.0);
  bfp3=c(0.0,2.00,0.0,1.79,2.00,0.0,0.0);
  rt=7.8;
  rl=24.0;
  LR_stat_vars<-list()
  #Setup variability coefficients (Section 31)
  if(propv$ctrl_state>0){

    if(propv$ctrl_state==4){
      #Mode of variability coefficients (Section 33)
      LR_stat_vars$kdv=propv$mode_variability;
      LR_stat_vars$ws = LR_stat_vars$kdv>=20;
      if(LR_stat_vars$ws){
        LR_stat_vars$kdv=LR_stat_vars$kdv-20;
      }
      LR_stat_vars$w1 = LR_stat_vars$kdv>=10;
      if(LR_stat_vars$w1){
        LR_stat_vars$kdv=LR_stat_vars$kdv - 10;
      }
      if(LR_stat_vars$kdv<0 | LR_stat_vars$kdv>3){
        LR_stat_vars$kdv=0;
        prop$error_id=max(prop$error_id,2);
        prop$errmsg = rbind(prop$errmsg, 'Warning: Mode of variablity out of range, using 0');
      }

      #Frequency Coefficients (Section 34)
      q=log(0.133*prop$wave_number);
      LR_stat_vars$gm=LR_stat_vars$cfm1+LR_stat_vars$cfm2/(((LR_stat_vars$cfm3*q)^2.0)+1.0);
      LR_stat_vars$gp=LR_stat_vars$cfp1+LR_stat_vars$cfp2/(((LR_stat_vars$cfp3*q)^2.0)+1.0);

      #System Coefficients (Section 35)
      LR_stat_vars$dexa=sqrt(18e6*prop$effective_height[1])+sqrt(18e6*prop$effective_height[2]) + ((575.7e12/prop$wave_number)^(1/3));
    }
    if(propv$ctrl_state==3){
      #Frequency Coefficients (Section 34)
      q=log(0.133*prop$wave_number);
      gm=cfm1+cfm2/(((cfm3*q)^2.0)+1.0);
      gp=cfp1+cfp2/(((cfp3*q)^2.0)+1.0);

      #System Coefficients (Section 35)
      dexa=sqrt(18e6*prop$effective_height[1])+sqrt(18e6*prop$effective_height[2]) + ((575.7e12/prop$wave_number)^(1/3));
    }
    if(propv$ctrl_state== 2){
      #System Coefficients (Section 35)
      dexa=sqrt(18e6*prop$effective_height[1])+sqrt(18e6*prop$effective_height[2]) + ((575.7e12/prop$wave_number)^(1/3));
    }

    if(propv$ctrl_state== 5){
      if(propv$climate_id<=0 | propv$climate_id>7){
        propv$climate_id = 5; #assume temperate

        prop$error_id=max(prop$error_id,2);
        prop$errmsg = rbind(prop$errmsg, 'Warning: Climate out of range, assuming Continental Temperate (climate code set to 5)');
      }
      #climatic coefficients (Section 32)
      LR_stat_vars$cv1 = bv1[propv$climate_id];
      LR_stat_vars$cv2 = bv2[propv$climate_id];
      LR_stat_vars$yv1 = xv1[propv$climate_id];
      LR_stat_vars$yv2 = xv2[propv$climate_id];
      LR_stat_vars$yv3 = xv3[propv$climate_id];
      LR_stat_vars$csm1=bsm1[propv$climate_id];
      LR_stat_vars$csm2=bsm2[propv$climate_id];
      LR_stat_vars$ysm1=xsm1[propv$climate_id];
      LR_stat_vars$ysm2=xsm2[propv$climate_id];
      LR_stat_vars$ysm3=xsm3[propv$climate_id];
      LR_stat_vars$csp1=bsp1[propv$climate_id];
      LR_stat_vars$csp2=bsp2[propv$climate_id];
      LR_stat_vars$ysp1=xsp1[propv$climate_id];
      LR_stat_vars$ysp2=xsp2[propv$climate_id];
      LR_stat_vars$ysp3=xsp3[propv$climate_id];
      LR_stat_vars$csd1=bsd1[propv$climate_id];
      LR_stat_vars$zd  =bzd1[propv$climate_id];
      LR_stat_vars$cfm1=bfm1[propv$climate_id];
      LR_stat_vars$cfm2=bfm2[propv$climate_id];
      LR_stat_vars$cfm3=bfm3[propv$climate_id];
      LR_stat_vars$cfp1=bfp1[propv$climate_id];
      LR_stat_vars$cfp2=bfp2[propv$climate_id];
      LR_stat_vars$cfp3=bfp3[propv$climate_id];

      #Mode of variability coefficients (Section 33)
      LR_stat_vars$kdv=propv$mode_variability;
      LR_stat_vars$ws = LR_stat_vars$kdv>=20;
      if(LR_stat_vars$ws){
        LR_stat_vars$kdv=LR_stat_vars$kdv-20;
      }
      LR_stat_vars$w1 = LR_stat_vars$kdv>=10;
      if(LR_stat_vars$w1){
        LR_stat_vars$kdv=LR_stat_vars$kdv-10;
      }
      if(LR_stat_vars$kdv<0 | LR_stat_vars$kdv>3){
        LR_stat_vars$kdv=0;
        prop$error_id=max(prop$error_id,2);
        prop$errmsg = rbind(prop$errmsg, 'Warning: Mode of variablity out of range, setting to 0\n');
      }

      #Frequency Coefficients (Section 34)
      q=log(0.133*prop$wave_number);
      LR_stat_vars$gm=LR_stat_vars$cfm1+LR_stat_vars$cfm2/(((LR_stat_vars$cfm3*q)^2)+1);
      LR_stat_vars$gp=LR_stat_vars$cfp1+LR_stat_vars$cfp2/(((LR_stat_vars$cfp3*q)^2)+1);

      #System Coefficients (Section 35)
      LR_stat_vars$dexa=sqrt(18e6*prop$effective_height[1])+sqrt(18e6*prop$effective_height[2]) + ((575.7e12/prop$wave_number)^(1/3));
    }
  }

  #Distance Coefficients (Section 36)
  if(prop$dist<LR_stat_vars$dexa){
    LR_stat_vars$de=130e3*prop$dist/LR_stat_vars$dexa;
  }else{
    LR_stat_vars$de=130e3+prop$dist-LR_stat_vars$dexa;
  }
  LR_stat_vars$vmd =  fCurve(LR_stat_vars$cv1, LR_stat_vars$cv2, LR_stat_vars$yv1, LR_stat_vars$yv2, LR_stat_vars$yv3, LR_stat_vars$de);
  LR_stat_vars$sgtm = fCurve(LR_stat_vars$csm1,LR_stat_vars$csm2,LR_stat_vars$ysm1,LR_stat_vars$ysm2,LR_stat_vars$ysm3,LR_stat_vars$de) * LR_stat_vars$gm;
  LR_stat_vars$sgtp = fCurve(LR_stat_vars$csp1,LR_stat_vars$csp2,LR_stat_vars$ysp1,LR_stat_vars$ysp2,LR_stat_vars$ysp3,LR_stat_vars$de) * LR_stat_vars$gp;
  LR_stat_vars$sgtd = LR_stat_vars$sgtp*LR_stat_vars$csd1;
  LR_stat_vars$tgtd = (LR_stat_vars$sgtp-LR_stat_vars$sgtd)*LR_stat_vars$zd;
  if(LR_stat_vars$w1){
    LR_stat_vars$sgl=0.0;
  }else{
    q=(1-0.8*exp(-prop$dist/50e3))*prop$terrain_irregularity*prop$wave_number;
    LR_stat_vars$sgl=10*q/(q+13);
  }
  if(LR_stat_vars$ws){
    LR_stat_vars$vs0=0;
  }else{
    LR_stat_vars$vs0=((5+3*exp(-LR_stat_vars$de/100e3))^2);
  }

  propv$ctrl_state=0;


  #correct normal deviations (Section 37)
  zt=zzt;
  zl=zzl;
  zc=zzc;

  if(LR_stat_vars$kdv== 0){
    zt=zc;
    zl=zc;
  }
  if(LR_stat_vars$kdv== 1){
    zl=zc;
  }
  if(LR_stat_vars$kdv==2){
    zl=zt;
  }

  if(abs(zt)>3.1 | abs(zl)>3.1 | abs(zc)>3.1){
    prop$error_id=max(prop$error_id,1);
    prop$errmsg = rbind(prop$errmsg, 'Warning: Normal deviates out of range for propagation loss model;\n');
  }


  #resolve standard deviations  (section 38)
  if(zt<0.0){
    sgt=LR_stat_vars$sgtm;
  }else if(zt<=LR_stat_vars$zd){
    sgt=LR_stat_vars$sgtp;
  }else{
    sgt=LR_stat_vars$sgtd+LR_stat_vars$tgtd/zt;
  }
  vs=LR_stat_vars$vs0+(sgt*zt)^2.0/(rt+zc^2)+(LR_stat_vars$sgl*zl)^2.0/(rl+zc^2);

  #Resolve deviations yr, yc. (Section 39)

  if(LR_stat_vars$kdv==0){
    yr=0.0;
    propv$sgc_confidence=sqrt(sgt^2+LR_stat_vars$sgl^2+vs);
  }
  if(LR_stat_vars$kdv==1){
    yr=sgt*zt;
    propv$sgc_confidence=sqrt(LR_stat_vars$sgl^2+vs);
  }
  if(LR_stat_vars$kdv==2){
    yr=sqrt(sgt^2+LR_stat_vars$sgl^2)*zt;
    propv$sgc_confidence=sqrt(vs);
  }
  if(!LR_stat_vars$kdv%in%c(0,1,2)){
    yr=sgt*zt+LR_stat_vars$sgl*zl;
    propv$sgc_confidence=sqrt(vs);
  }

  avarv=prop$reference_att-LR_stat_vars$vmd-yr-propv$sgc_confidence*zc;
  if(avarv<0.0){
    avarv=avarv*(29.0-avarv)/(29.0-10.0*avarv);
  }

  returnme<-list(avarv, prop,propv)
  names(returnme)<-c("avarv", "prop","propv")
  returnme
}


#' Area Mode Irregular Terrain Modeling
#'
#' Returns path loss in Longley Rice area mode.
#' @param  ModVar One of:
#'    \itemize{
#'   \item  0 - Single: pctConf is "Time/Situation/Location", pctTime, pctLoc not used
#'   \item     1 - Individual: pctTime is "Situation/Location", pctConf is "Confidence", pctLoc not used
#'   \item      2 - Mobile: pctTime is "Time/Locations (Reliability)", pctConf is "Confidence", pctLoc not used
#'   \item      3 - Broadcast: pctTime is "Time", pctLoc is "Location", pctConf is "Confidence"
#'   }
#' @param deltaH Terrain irregularity
#' @param tht_m Transmit antenna height above ground, m
#' @param rht_m Receive antenna height above ground, m
#' @param dist_km Distance to calculate db loss (radius dist in km from tower)
#' @param TSiteCriteria  0 - random, 1 - careful, 2 - very careful
#' @param RSiteCriteria  0 - random, 1 - careful, 2 - very careful
#' @param eps_dielect Soil dielectric
#' @param sgm_conductivity Surface conductivity
#' @param eno_ns_surfref Surface refractivity
#' @param frq_mhz Frequency to calculate loss at (Hz)
#' @param radio_climate 1-Equatorial, 2-Continental Subtropical, 3-Maritime Tropical, 4-Desert, 5-Continental Temperate, 6-Maritime Temperate, Over Land, 7-Maritime Temperate, Over Sea
#' @param pol Polarization. 0-Horizontal, 1-Vertical
#' @param pctTime Varies. (see parameter ModVar)
#' @param pctLoc Varies. (see parameter ModVar)
#' @param pctConf Varies. (see parameter ModVar)
#' @return Path loss (dB) and needed calculations for path loss
#' @examples
#'
#' ModVar = 3;
#' deltaH = 90;
#' tht_m = 100;
#' rht_m = 10;
#' dist_km = 20;
#' TSiteCriteria = 0;
#' RSiteCriteria = 0;
#' eps_dielect = 15;
#' sgm_conductivity = 0.005;
#' eno_ns_surfref = 301;
#' frq_mhz = 145;
#' radio_climate = 1;
#' pol = 1;      #1 = vert
#' pctTime = 0.5;
#' pctLoc = 0.5;
#' pctConf = 0.9;
#'
#' areaT(ModVar, deltaH, tht_m, rht_m, dist_km, TSiteCriteria, RSiteCriteria,
#' eps_dielect, sgm_conductivity, eno_ns_surfref,frq_mhz, radio_climate, pol, pctTime, pctLoc,
#' pctConf)$dbloss
#'
#' @export

areaT<-function(ModVar, deltaH, tht_m, rht_m, dist_km, TSiteCriteria,
                RSiteCriteria, eps_dielect, sgm_conductivity, eno_ns_surfref,frq_mhz,
                radio_climate, pol, pctTime, pctLoc, pctConf){
#tStart = tic; #?
# pol: 0-Horizontal, 1-Vertical
# TSiteCriteria, RSiteCriteria:
#		   0 - random, 1 - careful, 2 - very careful
# radio_climate: 1-Equatorial, 2-Continental Subtropical, 3-Maritime Tropical,
#                4-Desert, 5-Continental Temperate, 6-Maritime Temperate, Over Land,
#                7-Maritime Temperate, Over Sea
#
# pctTime, pctLoc, pctConf: .01 to .99
# errnum: 0- No Error.
#         1- Warning: Some parameters are nearly out of range.
#                     Results should be used with caution.
#         2- Note: Default parameters have been substituted for impossible ones.
#         3- Warning: A combination of parameters is out of range.
#                     Results are probably invalid.
#         Other-  Warning: Some parameters are out of range.
#                          Results are probably invalid.
# NOTE: strmode is not used at this time.

  strmode = NA;
  errnum = 0;

  list2env(genProps(),envir = environment())

  kst =  c(TSiteCriteria, RSiteCriteria)

  zt = qnorm(pctTime);
  zl = qnorm(pctLoc);
  zc = qnorm(pctConf);
  eps = eps_dielect;
  sgm = sgm_conductivity;
  eno = eno_ns_surfref;
  prop$terrain_irregularity = deltaH;

  prop$antenna_height = c(tht_m, rht_m);

  propv$climate_id =radio_climate;
  prop$surface_refractivity = eno;
  prop$error_id = 0;
  ivar = ModVar;
  ipol =  pol;

  prop=fLongleyRicePropSetup(frq_mhz, 0.0, eno, ipol, eps, sgm, prop);
  list2env(fLongleyRiceArea(kst, propv$climate_id, ivar, prop, propv),envir = environment())

  if(propv$ctrl_state<1) propv$ctrl_state = 1;

  list2env(fLongleyRiceProp(dist_km * 1000.0, prop, propa),envir = environment()) #[prop, propa]

  fs = 32.45 + 20.0 * log10(frq_mhz) + 20.0 * log10(prop$dist / 1000.0);

  templist<-fAttVariability(zt, zl, zc, prop, propv); #[xlb, prop,propv] =
  names(templist)<-c("xlb","prop","propv")
  list2env(templist,envir = environment())
  xlb = xlb + fs;
  dbloss = xlb;

  #if(prop$error_id==0){
  #  errnum = 0;
  #}else{
  #  errnum = prop$error_id;
  #}

  #Save common variables for quick recalculation
  #obj.struct_Internal.prop = prop;
  #obj.struct_Internal.propa = propa;
  #obj.struct_Internal.propv = propv;

  #toc(tStart)
  returnme<-list(dbloss,prop,propa,propv)
  names(returnme)<-c("dbloss","prop","propa","propv")
  returnme
}

#tStart = tic; #?
# pol: 0-Horizontal, 1-Vertical
# TSiteCriteria, RSiteCriteria:
#		   0 - random, 1 - careful, 2 - very careful
# radio_climate: 1-Equatorial, 2-Continental Subtropical, 3-Maritime Tropical,
#                4-Desert, 5-Continental Temperate, 6-Maritime Temperate, Over Land,
#                7-Maritime Temperate, Over Sea
# ModVar: 0 - Single: pctConf is "Time/Situation/Location", pctTime, pctLoc not used
#         1 - Individual: pctTime is "Situation/Location", pctConf is "Confidence", pctLoc not used
#         2 - Mobile: pctTime is "Time/Locations (Reliability)", pctConf is "Confidence", pctLoc not used
#         3 - Broadcast: pctTime is "Time", pctLoc is "Location", pctConf is "Confidence"
# pctTime, pctLoc, pctConf: .01 to .99
# errnum: 0- No Error.
#         1- Warning: Some parameters are nearly out of range.
#                     Results should be used with caution.
#         2- Note: Default parameters have been substituted for impossible ones.
#         3- Warning: A combination of parameters is out of range.
#                     Results are probably invalid.
#         Other-  Warning: Some parameters are out of range.
#                          Results are probably invalid.
# NOTE: strmode is not used at this time.







polyfit<-function (x, y, n = 1)
{
  if (!is.numeric(x) || !is.numeric(y))
    stop("Arguments x and y must be numeric.")
  if (length(x) != length(y))
    stop("Vectors/matrices x and y must be of same length.")
  if (is.null(n) || n < 0 || ceiling(n) != floor(n))
    stop("Degree n must be a non-negative integer.")
  x <- x[1:length(x)]
  y <- y[1:length(y)]
  A <- outer(x, seq(n, 0), "^")
  p <- qr.solve(A, y)
  return(p)
}


fLsqFit<-function(z, x1, x2){ #function [z0,zn] =
  #z=path_profile,x1=xl[1],x2=0.9*prop$horizon_dist[1]
  xstart = floor(x1/z[2]);
  xend = z[1] - floor(x2/z[2]);
  if (xend<0){
    xend = 0;
  }
  xend = z[1]-xend;
  xx = seq(xstart*z[2],xend*z[2],z[2])
  if (length(xx) == 1){
    z0 = z[3];
    zn = z0;
  }else{
    p = polyfit(xx,z[3+xstart:xend],1); #linear least square curve fit;
    #p = polyfit(xx,z(3:end),1); #linear least square curve fit;
    z0 = p[2];
    zn = p[1]*z[1]*z[2] + p[2];
  }
  returnme<-list(z0,zn)
  names(returnme)<-c("z0","zn")
  returnme
}



fHorizons<-function(path_profile, prop){ #prop
  end<-length(path_profile)
  np=path_profile[1]; #path_profile is the terrain profile, np = num points, xi = length of each increment
  xi=path_profile[2];
  za=path_profile[3]+prop$antenna_height[1];      #initial antenna height
  zb=path_profile[end]+prop$antenna_height[2];   #final antenna height
  qc=0.5*prop$earth_eff_curvature;
  q=qc*prop$dist;
  prop$horizon_elev_angle[2]=(zb-za)/prop$dist;
  prop$horizon_elev_angle[1]=prop$horizon_elev_angle[2]-q;
  prop$horizon_elev_angle[2]=-prop$horizon_elev_angle[2]-q;
  prop$horizon_dist[1]=prop$dist;
  prop$horizon_dist[2]=prop$dist;
  if(np>=2){
    sa=0.0;
    sb=prop$dist;
    wq=1;
    for (i in 2:(np)){
      sa= sa+xi;
      sb=sb-xi;
      q=path_profile[i+2]-(qc*sa+prop$horizon_elev_angle[1])*sa-za;
      if(q>0.0){
        prop$horizon_elev_angle[1]=prop$horizon_elev_angle[1]+q/sa;
        prop$horizon_dist[1]=sa;
        wq=0;
      }
      if(!wq){
        q=path_profile(i+2)-(qc*sb+prop$horizon_elev_angle[2])*sb-zb;
        if(q>0.0){
          prop$horizon_elev_angle[2]=prop$horizon_elev_angle[2]+q/sb;
          prop$horizon_dist[2]=sb;
        }
      }
    }
  }
  return(prop)
}


fLongleyRicePointtoPoint<-function(path_profile, klimx, mdvarx,prop, propa, propv ){ #function [prop, propa, propv] =
  #path_profile=elev,klimx=propv$climate_id,mdvarx=propv$mode_variability,prop,propa,propv
  prop$dist=path_profile[1]*path_profile[2];
  prop = fHorizons(path_profile,prop);
  end<-length(path_profile)
  xl = c(0 ,0);
  for (j in 1:2){
    xl[j]=min(15.0*prop$antenna_height[j],0.1*prop$horizon_dist[j]);
  }
  xl[2]=prop$dist-xl[2];
  prop$terrain_irregularity=fCalcTerrainIrregularity(path_profile,xl[1],xl[2]);

  if((prop$horizon_dist[1]+prop$horizon_dist[2])>1.5*prop$dist){
    templist<-fLsqFit(path_profile,xl[1],xl[2])
    names(templist)<-c("za","zb")
    list2env(templist,envir = environment()); #[za,zb]=

    prop$effective_height[1]=prop$antenna_height[1]+max(path_profile[3]-za,0);
    prop$effective_height[2]=prop$antenna_height[2]+max(path_profile[end]-zb,0);
    for (j in 1:2){
      prop$horizon_dist[j]=sqrt(2.0*prop$effective_height[j]/prop$earth_eff_curvature) *exp(-0.07*sqrt(prop$terrain_irregularity/max(prop$effective_height[j],5.0)));
    }
    q=prop$horizon_dist[1]+prop$horizon_dist[2];

    if(q<=prop$dist){
      q=((prop$dist/q)^2);
      for (j in 1:2){
        prop$effective_height[j]=prop$effective_height[j]*q;
        prop$horizon_dist[j]=sqrt(2.0*prop$effective_height[j]/prop$earth_eff_curvature) *exp(-0.07*sqrt(prop$terrain_irregularity/max(prop$effective_height[j],5.0)));
      }
    }
    for (j in 1:2){
      q=sqrt(2.0*prop$effective_height[j]/prop$earth_eff_curvature);
      prop$horizon_elev_angle[j] = (0.65*prop$terrain_irregularity*(q/prop$horizon_dist[j]-1.0)-2.0 *  prop$effective_height[j])/q;
    }

  }else{
    za=fLsqFit(path_profile,xl[1],0.9*prop$horizon_dist[1])[1]; #[za,~] =
    zb=fLsqFit(path_profile,prop$dist-0.9*prop$horizon_dist[2],xl[2])[2]; #[~,zb] =
    prop$effective_height[1] = prop$antenna_height[1]+max(path_profile[3]-za,0);
    prop$effective_height[2] = prop$antenna_height[2]+max(path_profile[end]-zb,0);
  }
  prop$mode = -1; #point to point
  propv$ctrl_state = max(propv$ctrl_state,3);
  if(mdvarx>= 0){
    propv$mode_variability = mdvarx;
    propv$ctrl_state = max(propv$ctrl_state,4);
  }
  if(klimx>0){
    propv$climate_id = klimx;
    propv$ctrl_state = 5;
  }
  list2env(fLongleyRiceProp(0,prop,propa),envir = environment()); #[prop, propa] =
  returnme<-(list(prop,propa,propv))
  names(returnme)<-c("prop","propa","propv")
  returnme
}

#' Point to Point ITM
#'
#' Returns path loss in Longley Rice point-to-point mode. Best for rural areas.
#' @param struct_Input Named list object of input parameters:
#' \itemize{
#'   \item Frequency - Frequency to calculate loss at (Hz)
#'   \item Elevation - terrain elevation profile, (list of points) (m)
#'   \item Resolution - terrain input resolution (distance b/t points) (m)
#'   \item TX_Height - Transmit antenna height above ground (m)
#'   \item RX_Height - Recieve antenna height above ground (m)
#'   \item eps - Soil dielectric
#'   \item sgm - Surface conductivity
#'   \item surfref - Surface refractivity
#'   \item Climate -  Climate, 1-Equitorial, 2-Continental Subtropical, 3-Maritime Tropical, 4-Desert
#                            5-Continental Temperate, 6-Maritime Temperate Over Lane 7-Maritime Temperate Over Sea
#'   \item Polarization - 1 is vertical, 0 is horizontal.
#'   \item Confidence - confidence for statistical analysis (.01 to .99)
#'   \item Reliability - Reliability to calculate statistics for (.01 to .99)
#' }

#' @return Path loss (dB), error ID, error message, and mode.
#' @examples
#'
#' #commented below is an example of how to get an elevation profile in R.
#' #library(elevatr)
#' #library(sp)
#' #set.seed(65.7)
#' #examp_df <- data.frame(x = runif(3, min = -73, max = -72.5), y = runif(3, min = 42,  max = 43))
#' #prj_dd <- "+init=EPSG:4326"
#' #cats <- data.frame(category = c("H", "M", "L"))
#' #examp_df2 <- data.frame(examp_df, cats)
#' #examp_sp <- SpatialPoints(examp_df, proj4string = CRS(prj_dd))
#' #examp_spdf <- SpatialPointsDataFrame(examp_sp, data = cats)
#' #df_elev_epqs <- get_elev_point(examp_df, prj = prj_dd, src = "epqs")
#' #Elevation<-df_elev_epqs$elevation
#'
#' #These are the values returned above:
#' Elevation<-c(207.81, 198.95, 306.15)
#'
#' #Build the input list
#' struct_Input<-list()
#' struct_Input$Frequency<-120*1000000 #Frequency to calculate loss at (Hz)
#' struct_Input$Elevation<-Elevation #terrain elevation profile, (list of points) (m)
#' struct_Input$Resolution<-40000 #terrain input resolution (distance b/t points) (m)
#' struct_Input$TX_Height<-3 #Transmit antenna height above ground (m)
#' struct_Input$RX_Height<-100 #Recieve antenna height above ground (m)
#' struct_Input$eps<-15 #Soil dielectric
#' struct_Input$sgm<-.005 #Surface conductivity
#' struct_Input$surfref<-301 #Surface refractivity
#' struct_Input$Climate<-5
#' struct_Input$Polarization<-1 #1 is vertical, 0 is horizontal
#' struct_Input$Confidence<-.95 #confidence for statistical analysis
#' struct_Input$Reliability<-.95 #Reliability to calculate statistics for (.01 to .99)
#'
#' point_to_point(struct_Input)
#' @export


#Point to Point Longley-Rice
point_to_point<-function(struct_Input){      #function struct_Output         =
  #===================================================
  #Function Name:
  #    point_to_point
  #
  #Inputs:
  #    struct_Input - Input containing the fields:
  #       .Frequency - Frequency to calculate loss at (Hz)
  #       .Elevation - terrain elevation profile, (list of points) (m)
  #       .Resolution - terrain input resolution (distance b/t points) (m)
  #       .TX_Height - Transmit antenna height above ground (m)
  #       .RX_Height - Recieve antenna height above ground (m)
  #       .eps - Soil dielectric
  #       .sgm - Surface conductivity
  #       .surfref - Surface refractivity
  #       .Climate -  Climate, 1-Equitorial, 2-Continental Subtropical, 3-Maritime Tropical, 4-Desert
  #                            5-Continental Temperate, 6-Maritime Temperate Over Lane 7-Maritime Temperate Over Sea
  #       .Polarization - 1 is vertical, 0 is horizontal.
  #       .Confidence - confidence for statistical analysis (.01 to .99)
  #       .Reliability - Reliability to calculate statistics for (.01 to .99)
  #
  #Outputs:
  #    struct_Output - Output containing the fields:
  #       .dbloss - Final Output loss (dB)
  #       .error_ID - Error ID, 0 = no error, 1 = may be out of range, 3 = out of range
  #       .error_msg - Detailed error info (string)
  #
  #Functions Called:
  #    fLongleyRicePropSetup
  #    fLongleyRicePointtoPoint
  #    fAttVariability
  #
  #Description:
  #    Calculate path loss using the Longley-Rice Irregular terrain
  #    model.
  #===================================================



  frq_mhz = struct_Input$Frequency /1e6;               #input should be in Hz
  elev = struct_Input$Elevation;                       #terrain input in meters
  Resolution = struct_Input$Resolution;                #terrain input resolution
  tht_m = struct_Input$TX_Height;                      #Transmit antenna height above ground
  rht_m = struct_Input$RX_Height;                      #Recieve antenna height above ground
  eps_dielect = struct_Input$eps;                      #Soil dielectric
  sgm_conductivity = struct_Input$sgm;
  eno_ns_surfref = struct_Input$surfref;

  radio_climate = struct_Input$Climate;                #Climate (1-7)
  pol = struct_Input$Polarization;                     #1 is vertical, 0 is horizontal.

  conf = struct_Input$Confidence;                      #confidence for statistical analysis
  rel = struct_Input$Reliability;                      #link reliability for statistical analysis



  elev = c(length(elev)-1, Resolution, elev)

  #[prop, propa, propv] = obj.genProps();
  #Initialize Propagation structures for shared variables
  prop = list();
  propa = list();
  propv = list();
  prop$errmsg = list(); #Initialize Error capture.

  zsys=0;
  strmode = list();
  prop$antenna_height = c(tht_m, rht_m)

  propv$climate_id = radio_climate;
  prop$error_id = 0;
  propv$ctrl_state = 5;
  prop$mode = -1; #point to point
  zc = qnorm(conf);
  zr = qnorm(rel);
  np = elev[1];
  dkm = (elev[2] * elev[1]) / 1000.0;
  xkm = elev[2] / 1000.0;
  eno = eno_ns_surfref;
  enso = 0.0;
  q = enso;
  if(q<=0.0){
    if (np > 3){
      ja = 3.0 ;#+ 0.1 * elev[1];
      jb = np;#np - ja + 6;
      for (i in (ja):(jb)){         #mod from (ja-1):(jb-1) ???
        zsys = zsys + elev[floor(i)];
      }

      zsys =  zsys/(jb-ja+1);
    }else{
      zsys = mean(elev[3:length(elev)]);

    }

    q = eno;
  }
  propv$mode_variability = 12; #Climate mode of variability, used for area mode when dealing with location, time, situation variability
  prop = fLongleyRicePropSetup(frq_mhz,zsys,q,pol,eps_dielect,sgm_conductivity,prop);
  list2env(fLongleyRicePointtoPoint(elev,propv$climate_id,propv$mode_variability,prop,propa,propv),envir = environment()); #[prop, propa, propv]  =
  fs = 32.45 + 20.0 * log10(frq_mhz) + 20.0 * log10(prop$dist / 1000.0);


  q = prop$dist - propa$total_horizon_dist;
  if(q<0.0){
    strmode='Line-Of-Sight Mode';
  }else{
    if(floor(q)==0.0){
      strmode = paste0(strmode, 'Single Horizon')

    }else if(floor(q)>0.0){
      strmode = paste0(strmode, 'Horizon')

    }
    if(prop$dist<=propa$los_distance | prop$dist <= propa$scatter_dist){
      strmode = paste0(strmode, ', Diffraction Dominant')

    }else if(prop$dist>propa$scatter_dist){
      strmode = paste0(strmode, ', Troposcatter Dominant')

    }
  }

  templist<-fAttVariability(zr,0.0,zc,prop,propv); #[dbloss, prop, propv] =
  names(templist)<-c("dbloss", "prop", "propv")
  list2env(templist,envir= environment())
  struct_Output<-list()
  struct_Output$dbloss = dbloss + fs;
  #   fprintf([num2str(struct_Output$dbloss) '\n' num2str(struct_Output$dbloss) '\n' num2str(struct_Output$dbloss) '\n' num2str(struct_Output$dbloss) '\n' num2str(struct_Output$dbloss) '\n' num2str(struct_Output$dbloss) '\n'  ])
  struct_Output$error_ID = prop$error_id;
  struct_Output$error_msg = prop$errmsg;
  struct_Output$mode = strmode;

  #   if prop$error_id~=0;   fprintf(prop$errmsg); end;
  #   fprintf([strmode '\n']);

  struct_Output
}
