#'  Rounde
#'
#'  Rounding function - SPlus rounding function uses the nearest even number rule
#' @keywords zscore igrowup who
#' @details WHO Child Growth Standards Department of Nutrition for Health and Development Last modified on 07/10/2013-Developed using R version 3.0.1
#' @note This code conrcerns the standard approach for the prevalences, i.e. the calculation of the prevalences takes into account all the valid (non-missing) z-scores for each of the indicators.
#' @export
#' @examples
#' x<-c(1,1)
#' rounde(x)

rounde <- function(x,digits=0) {
  expo<-10^digits
  return(ifelse(abs(x*expo) - floor(abs(x*expo)) < 0.5, sign(x*expo) * floor(abs(x*expo)), sign(x*expo) * (floor(abs(x*expo)) + 1))/expo)
}

#'  Prevph.L
#'
#'  Z-score prevalence calc Without oedema (for head circumference-for-age and arm circumferecen-for-age)
#' @keywords zscore igrowup who
#' @details WHO Child Growth Standards Department of Nutrition for Health and Development Last modified on 07/10/2013-Developed using R version 3.0.1
#' @note This code conrcerns the standard approach for the prevalences, i.e. the calculation of the prevalences takes into account all the valid (non-missing) z-scores for each of the indicators.
#' @export
#' @examples
#' a<-c(1,1)
#' x<-c(1,1)
#' w<-c(1,1)
#' prevph.L(a,x,w)

prevph.L <- function(a,x,w) {
  ph <- sum((x > a)*w,na.rm=T)/sum((!is.na(x))*w,na.rm=T)
  aux <- 1.96*sqrt(ph*(1-ph)/sum((!is.na(x))*w,na.rm=T))+(1/(2*sum((!is.na(x))*w,na.rm=T)))
  vec <- c(rounde(sum((!is.na(x))*w,na.rm=T),digits=0),rounde(c(ph, max(0,ph-aux), ph+aux)*100,digits=1))
  return(vec)
}

#'  Prevph
#'
#'  Z-score prevalence calc With oedema (only for weight-for-length/height and bmi-for-age)
#' @keywords zscore igrowup who
#' @details WHO Child Growth Standards Department of Nutrition for Health and Development Last modified on 07/10/2013-Developed using R version 3.0.1
#' @note This code conrcerns the standard approach for the prevalences, i.e. the calculation of the prevalences takes into account all the valid (non-missing) z-scores for each of the indicators.
#' @export
#' @examples
#' a<-c(1,1)
#' x<-c(1,1)
#' w<-c(1,1)
#' f<-"y"
#' prevph(a,x,w,f)

prevph <- function(a,x,w,f) {
  f<-as.character(f)
  ph <- sum((x > a)*w,na.rm=T)/sum(((!is.na(x)) | (f=="y"))*w,na.rm=T)
  aux <- 1.96*sqrt(ph*(1-ph)/sum(((!is.na(x)) | (f=="y"))*w,na.rm=T))+(1/(2*sum(((!is.na(x)) | (f=="y"))*w,na.rm=T)))
  vec <- c(rounde(sum(((!is.na(x)) | (f=="y"))*w,na.rm=T),digits=0),rounde(c(ph, max(0,ph-aux), ph+aux)*100,digits=1))
  return(vec)
}

#'  Prevnh.L
#'
#'  Prevalence calculation for the lower bound and corresponding 95% C.I. calc Without oedema (for length/height-for-age)
#' @keywords zscore igrowup who
#' @details WHO Child Growth Standards Department of Nutrition for Health and Development Last modified on 07/10/2013-Developed using R version 3.0.1
#' @note This code conrcerns the standard approach for the prevalences, i.e. the calculation of the prevalences takes into account all the valid (non-missing) z-scores for each of the indicators.
#' @export
#' @examples
#' a<-c(1,1)
#' x<-c(1,1)
#' w<-c(1,1)
#' prevnh.L(a,x,w)

prevnh.L <- function(a,x,w) {
  ph <- sum((x < a)*w,na.rm=T)/sum((!is.na(x))*w,na.rm=T)
  aux <- 1.96*sqrt(ph*(1-ph)/sum((!is.na(x))*w,na.rm=T))+(1/(2*sum((!is.na(x))*w,na.rm=T)))
  vec <- c(rounde(sum((!is.na(x))*w,na.rm=T),digits=0),rounde(c(ph, max(0,ph-aux), ph+aux)*100,digits=1))
  return(vec)
}

#'  prevnh
#'
#'  Prevalence calculation for the lower bound and corresponding 95% C.I. calc With oedema  (for all weight-related indicators)
#' @keywords zscore igrowup who
#' @details WHO Child Growth Standards Department of Nutrition for Health and Development Last modified on 07/10/2013-Developed using R version 3.0.1
#' @note This code conrcerns the standard approach for the prevalences, i.e. the calculation of the prevalences takes into account all the valid (non-missing) z-scores for each of the indicators.
#' @export
#' @examples
#' a<-c(1,1)
#' x<-c(1,1)
#' w<-c(1,1)
#' f<-"y"
#' prenph(a,x,w,f)

prevnh <- function(a,x,w,f) {
  f<-as.character(f)
  ph <- sum((x < a | f=="y")*w,na.rm=T)/sum(((!is.na(x)) | (f=="y"))*w,na.rm=T)
  aux <- 1.96*sqrt(ph*(1-ph)/sum(((!is.na(x)) | (f=="y"))*w,na.rm=T))+(1/(2*sum(((!is.na(x)) | (f=="y"))*w,na.rm=T)))
  vec <- c(rounde(sum(((!is.na(x)) | (f=="y"))*w,na.rm=T),digits=0),rounde(c(ph, max(0,ph-aux), ph+aux)*100,digits=1))
  return(vec)
}


#'  wmean
#'
#'  Weighted mean
#' @keywords zscore igrowup who
#' @details WHO Child Growth Standards Department of Nutrition for Health and Development Last modified on 07/10/2013-Developed using R version 3.0.1
#' @note This code conrcerns the standard approach for the prevalences, i.e. the calculation of the prevalences takes into account all the valid (non-missing) z-scores for each of the indicators.
#' @export
#' @examples
#' x<-c(1,1)
#' w<-c(1,1)
#' wmean(x,w)

wmean <- function(x,w) {
  return(rounde(sum(x*w,na.rm=T)/sum(w[!is.na(x)]),digits=2) )
  }

#'  wsd
#'
#'  Weighted standard deviation
#' @keywords zscore igrowup who
#' @details WHO Child Growth Standards Department of Nutrition for Health and Development Last modified on 07/10/2013-Developed using R version 3.0.1
#' @note This code conrcerns the standard approach for the prevalences, i.e. the calculation of the prevalences takes into account all the valid (non-missing) z-scores for each of the indicators.
#' @export
#' @examples
#' x<-c(1,1)
#' w<-c(1,1)
#' wsd(x,w)

wsd <- function(x,w) {
  mh <- sum(x*w,na.rm=T)/sum((!is.na(x))*w,na.rm=T)
  sdh<-ifelse(length(x[!is.na(x)])>0,rounde(sqrt(sum(((x-mh)^2)*w,na.rm=T)/(sum((!is.na(x))*w,na.rm=T) - 1)),digits=2),NA)
  return( sdh )
}


#'  calc.zlen
#'
#'  Function for calculating individual Length-for-age z-scores
#' @keywords zscore igrowup who
#' @details WHO Child Growth Standards Department of Nutrition for Health and Development Last modified on 07/10/2013-Developed using R version 3.0.1
#' @note This code conrcerns the standard approach for the prevalences, i.e. the calculation of the prevalences takes into account all the valid (non-missing) z-scores for each of the indicators.
#' @export
#' @examples
#' data("igrowupstd_AnXprs")
#' mat<-data.frame(age.days=c(1854,1854),sex=c(1,2),clenhei=c(110.4611,109.6604))
#' calc.zlen(mat,lenanthro)

calc.zlen<-function(mat,lenanthro){
  for(i in 1:length(mat$age.days)) {
    if(!is.na(mat$age.days[i]) & mat$age.days[i]>=0 & mat$age.days[i]<=1856) {
      l.val<-lenanthro$l[lenanthro$age==mat$age.days[i] & lenanthro$sex==mat$sex[i]]
      m.val<-lenanthro$m[lenanthro$age==mat$age.days[i] & lenanthro$sex==mat$sex[i]]
      s.val<-lenanthro$s[lenanthro$age==mat$age.days[i] & lenanthro$sex==mat$sex[i]]
      mat$zlen[i]<-(((mat$clenhei[i]/m.val)^l.val)-1)/(s.val*l.val)
    }	else mat$zlen[i]<- NA
  }
  return(mat)
}

#'  calc.zhc
#'
#'  Function for calculating individual Head circumference-for-age z-scores
#' @keywords zscore igrowup who
#' @details WHO Child Growth Standards Department of Nutrition for Health and Development Last modified on 07/10/2013-Developed using R version 3.0.1
#' @note This code conrcerns the standard approach for the prevalences, i.e. the calculation of the prevalences takes into account all the valid (non-missing) z-scores for each of the indicators.
#' @export

calc.zhc<-function(mat,hcanthro){

  for(i in 1:length(mat$age.days)) {

    if(!is.na(mat$age.days[i]) & mat$age.days[i]>=0 & mat$age.days[i]<=1856) {

      l.val<-hcanthro$l[hcanthro$age==mat$age.days[i] & hcanthro$sex==mat$sex[i]]
      m.val<-hcanthro$m[hcanthro$age==mat$age.days[i] & hcanthro$sex==mat$sex[i]]
      s.val<-hcanthro$s[hcanthro$age==mat$age.days[i] & hcanthro$sex==mat$sex[i]]
      mat$zhc[i]<-(((mat$headc[i]/m.val)^l.val)-1)/(s.val*l.val)

    }	else mat$zhc[i]<- NA

  }
  return(mat)
}

#'  calc.zwei
#'
#'  Function for calculating individual Weight-for-age z-scores
#' @keywords zscore igrowup who
#' @details WHO Child Growth Standards Department of Nutrition for Health and Development Last modified on 07/10/2013-Developed using R version 3.0.1
#' @note This code conrcerns the standard approach for the prevalences, i.e. the calculation of the prevalences takes into account all the valid (non-missing) z-scores for each of the indicators.
#' @export
calc.zwei<-function(mat,weianthro){

  for(i in 1:length(mat$age.days)) {

    if(!is.na(mat$age.days[i]) & mat$age.days[i]>=0 & mat$age.days[i]<=1856 & mat$oedema[i]!="y") {

      l.val<-weianthro$l[weianthro$age==mat$age.days[i] & weianthro$sex==mat$sex[i]]
      m.val<-weianthro$m[weianthro$age==mat$age.days[i] & weianthro$sex==mat$sex[i]]
      s.val<-weianthro$s[weianthro$age==mat$age.days[i] & weianthro$sex==mat$sex[i]]

      mat$zwei[i]<-(((mat$weight[i]/m.val)^l.val)-1)/(s.val*l.val)
      if(!is.na(mat$zwei[i]) & mat$zwei[i]>3) {
        sd3pos<- m.val*((1+l.val*s.val*3)^(1/l.val))
        sd23pos<- sd3pos- m.val*((1+l.val*s.val*2)^(1/l.val))
        mat$zwei[i]<- 3+((mat$weight[i]-sd3pos)/sd23pos)
      }
      if(!is.na(mat$zwei[i]) & mat$zwei[i]< (-3)) {
        sd3neg<- m.val*((1+l.val*s.val*(-3))**(1/l.val))
        sd23neg<- m.val*((1+l.val*s.val*(-2))**(1/l.val))-sd3neg
        mat$zwei[i]<- (-3)+((mat$weight[i]-sd3neg)/sd23neg)
      }

    } else mat$zwei[i]<-NA
  }
  return(mat)
}

#'  calc.zac
#'
#'  Function for calculating individual Arm circumference-for-age z-scores
#' @keywords zscore igrowup who
#' @details WHO Child Growth Standards Department of Nutrition for Health and Development Last modified on 07/10/2013-Developed using R version 3.0.1
#' @note This code conrcerns the standard approach for the prevalences, i.e. the calculation of the prevalences takes into account all the valid (non-missing) z-scores for each of the indicators.
#' @export

calc.zac<-function(mat,acanthro){

  for(i in 1:length(mat$age.days)) {

    if(!is.na(mat$age.days[i]) & mat$age.days[i]>=91 & mat$age.days[i]<=1856) {

      l.val<-acanthro$l[acanthro$age==mat$age.days[i] & acanthro$sex==mat$sex[i]]
      m.val<-acanthro$m[acanthro$age==mat$age.days[i] & acanthro$sex==mat$sex[i]]
      s.val<-acanthro$s[acanthro$age==mat$age.days[i] & acanthro$sex==mat$sex[i]]
      mat$zac[i]<-(((mat$armc[i]/m.val)^l.val)-1)/(s.val*l.val)
      if(!is.na(mat$zac[i]) & mat$zac[i]>3) {
        sd3pos<- m.val*((1+l.val*s.val*3)^(1/l.val))
        sd23pos<- sd3pos- m.val*((1+l.val*s.val*2)^(1/l.val))
        mat$zac[i]<- 3+((mat$armc[i]-sd3pos)/sd23pos)
      }
      if(!is.na(mat$zac[i]) & mat$zac[i]< (-3)) {
        sd3neg<- m.val*((1+l.val*s.val*(-3))**(1/l.val))
        sd23neg<- m.val*((1+l.val*s.val*(-2))**(1/l.val))-sd3neg
        mat$zac[i]<- (-3)+((mat$armc[i]-sd3neg)/sd23neg)
      }

    } else mat$zac[i]<-NA

  }
  return(mat)
}

#'  calc.zts
#'
#'  Function for calculating individual Triceps skinfold-for-age z-scores
#' @keywords zscore igrowup who
#' @details WHO Child Growth Standards Department of Nutrition for Health and Development Last modified on 07/10/2013-Developed using R version 3.0.1
#' @note This code conrcerns the standard approach for the prevalences, i.e. the calculation of the prevalences takes into account all the valid (non-missing) z-scores for each of the indicators.
#' @export

calc.zts<-function(mat,tsanthro){

  for(i in 1:length(mat$age.days)) {

    if(!is.na(mat$age.days[i]) & mat$age.days[i]>=91 & mat$age.days[i]<=1856) {

      l.val<-tsanthro$l[tsanthro$age==mat$age.days[i] & tsanthro$sex==mat$sex[i]]
      m.val<-tsanthro$m[tsanthro$age==mat$age.days[i] & tsanthro$sex==mat$sex[i]]
      s.val<-tsanthro$s[tsanthro$age==mat$age.days[i] & tsanthro$sex==mat$sex[i]]

      mat$zts[i]<-(((mat$triskin[i]/m.val)^l.val)-1)/(s.val*l.val)
      if(!is.na(mat$zts[i]) & mat$zts[i]>3) {
        sd3pos<- m.val*((1+l.val*s.val*3)^(1/l.val))
        sd23pos<- sd3pos- m.val*((1+l.val*s.val*2)^(1/l.val))
        mat$zts[i]<- 3+((mat$triskin[i]-sd3pos)/sd23pos)
      }
      if(!is.na(mat$zts[i]) & mat$zts[i]< (-3)) {
        sd3neg<- m.val*((1+l.val*s.val*(-3))**(1/l.val))
        sd23neg<- m.val*((1+l.val*s.val*(-2))**(1/l.val))-sd3neg
        mat$zts[i]<- (-3)+((mat$triskin[i]-sd3neg)/sd23neg)
      }

    } else mat$zts[i]<-NA

  }

  return(mat)
}

#'  calc.zss
#'
#'  Function for calculating individual Subscapular skinfold-for-age z-scores
#' @keywords zscore igrowup who
#' @details WHO Child Growth Standards Department of Nutrition for Health and Development Last modified on 07/10/2013-Developed using R version 3.0.1
#' @note This code conrcerns the standard approach for the prevalences, i.e. the calculation of the prevalences takes into account all the valid (non-missing) z-scores for each of the indicators.
#' @export

calc.zss<-function(mat,ssanthro){

  for(i in 1:length(mat$age.days)) {

    if(!is.na(mat$age.days[i]) & mat$age.days[i]>=91 & mat$age.days[i]<=1856) {

      l.val<-ssanthro$l[ssanthro$age==mat$age.days[i] & ssanthro$sex==mat$sex[i]]
      m.val<-ssanthro$m[ssanthro$age==mat$age.days[i] & ssanthro$sex==mat$sex[i]]
      s.val<-ssanthro$s[ssanthro$age==mat$age.days[i] & ssanthro$sex==mat$sex[i]]

      mat$zss[i]<-(((mat$subskin[i]/m.val)^l.val)-1)/(s.val*l.val)
      if(!is.na(mat$zss[i]) & mat$zss[i]>3) {
        sd3pos<- m.val*((1+l.val*s.val*3)^(1/l.val))
        sd23pos<- sd3pos- m.val*((1+l.val*s.val*2)^(1/l.val))
        mat$zss[i]<- 3+((mat$subskin[i]-sd3pos)/sd23pos)
      }
      if(!is.na(mat$zss[i]) & mat$zss[i]< (-3)) {
        sd3neg<- m.val*((1+l.val*s.val*(-3))**(1/l.val))
        sd23neg<- m.val*((1+l.val*s.val*(-2))**(1/l.val))-sd3neg
        mat$zss[i]<- (-3)+((mat$subskin[i]-sd3neg)/sd23neg)
      }

    } else mat$zss[i]<-NA
  }

  return(mat)
}

#'  calc.zwfl
#'
#'  Function for calculating individual Weight-for-length/height z-scores
#' @keywords zscore igrowup who
#' @details WHO Child Growth Standards Department of Nutrition for Health and Development Last modified on 07/10/2013-Developed using R version 3.0.1
#' @note This code conrcerns the standard approach for the prevalences, i.e. the calculation of the prevalences takes into account all the valid (non-missing) z-scores for each of the indicators.
#' @export

calc.zwfl<-function(mat,wflanthro,wfhanthro){

  for(i in 1:length(mat$age.days)) {

    mat$zwfl[i]<-NA

    if(mat$oedema[i]!="y") {

      if( (!is.na(mat$age.days[i]) & mat$age.days[i]<731) | (is.na(mat$age.days[i]) & !is.na(mat$l.h[i]) & (mat$l.h[i]=="l" | mat$l.h[i]=="L")) | (is.na(mat$age.days[i]) & is.na(mat$l.h[i]) & !is.na(mat$clenhei[i]) & mat$clenhei[i]<87) ) {

        if(!is.na(mat$clenhei[i]) & mat$clenhei[i]>=45 & mat$clenhei[i]<=110) {

          ### Interpolated l,m,s values

          low.len<-trunc(mat$clenhei[i]*10)/10
          upp.len<-trunc(mat$clenhei[i]*10+1)/10
          diff.len<-(mat$clenhei[i]-low.len)/0.1

          if(diff.len>0) {
            l.val<-wflanthro$l[wflanthro$length==low.len & wflanthro$sex==mat$sex[i]]+diff.len*( wflanthro$l[wflanthro$length==upp.len & wflanthro$sex==mat$sex[i]]-wflanthro$l[wflanthro$length==low.len & wflanthro$sex==mat$sex[i]] )
            m.val<-wflanthro$m[wflanthro$length==low.len & wflanthro$sex==mat$sex[i]]+diff.len*( wflanthro$m[wflanthro$length==upp.len & wflanthro$sex==mat$sex[i]]-wflanthro$m[wflanthro$length==low.len & wflanthro$sex==mat$sex[i]] )
            s.val<-wflanthro$s[wflanthro$length==low.len & wflanthro$sex==mat$sex[i]]+diff.len*( wflanthro$s[wflanthro$length==upp.len & wflanthro$sex==mat$sex[i]]-wflanthro$s[wflanthro$length==low.len & wflanthro$sex==mat$sex[i]] )
          } else {
            l.val<-wflanthro$l[wflanthro$length==low.len & wflanthro$sex==mat$sex[i]]
            m.val<-wflanthro$m[wflanthro$length==low.len & wflanthro$sex==mat$sex[i]]
            s.val<-wflanthro$s[wflanthro$length==low.len & wflanthro$sex==mat$sex[i]]
          }

          mat$zwfl[i]<-(((mat$weight[i]/m.val)^l.val)-1)/(s.val*l.val)
          if(!is.na(mat$zwfl[i]) & mat$zwfl[i]>3) {
            sd3pos<- m.val*((1+l.val*s.val*3)^(1/l.val))
            sd23pos<- sd3pos- m.val*((1+l.val*s.val*2)^(1/l.val))
            mat$zwfl[i]<- 3+((mat$weight[i]-sd3pos)/sd23pos)
          }
          if(!is.na(mat$zwfl[i]) & mat$zwfl[i]<(-3)) {
            sd3neg<- m.val*((1+l.val*s.val*(-3))**(1/l.val))
            sd23neg<- m.val*((1+l.val*s.val*(-2))**(1/l.val))-sd3neg
            mat$zwfl[i]<- (-3)-((sd3neg-mat$weight[i])/sd23neg)
          }
        }
      }

      else 		if( (!is.na(mat$age.days[i]) & mat$age.days[i]>=731) | (is.na(mat$age.days[i]) & !is.na(mat$l.h[i]) & (mat$l.h[i]=="h" | mat$l.h[i]=="H"))  | (is.na(mat$age.days[i]) & is.na(mat$l.h[i]) & !is.na(mat$clenhei[i]) & mat$clenhei[i]>=87) ) {

        if(!is.na(mat$clenhei[i]) & mat$clenhei[i]>=65 & mat$clenhei[i]<=120) {

          ### Interpolated l,m,s values

          low.len<-trunc(mat$clenhei[i]*10)/10
          upp.len<-trunc(mat$clenhei[i]*10+1)/10
          diff.len<-(mat$clenhei[i]-low.len)/0.1

          if(diff.len>0) {
            l.val<-wfhanthro$l[wfhanthro$height==low.len & wfhanthro$sex==mat$sex[i]]+diff.len*( wfhanthro$l[wfhanthro$height==upp.len & wfhanthro$sex==mat$sex[i]]-wfhanthro$l[wfhanthro$height==low.len & wfhanthro$sex==mat$sex[i]] )
            m.val<-wfhanthro$m[wfhanthro$height==low.len & wfhanthro$sex==mat$sex[i]]+diff.len*( wfhanthro$m[wfhanthro$height==upp.len & wfhanthro$sex==mat$sex[i]]-wfhanthro$m[wfhanthro$height==low.len & wfhanthro$sex==mat$sex[i]] )
            s.val<-wfhanthro$s[wfhanthro$height==low.len & wfhanthro$sex==mat$sex[i]]+diff.len*( wfhanthro$s[wfhanthro$height==upp.len & wfhanthro$sex==mat$sex[i]]-wfhanthro$s[wfhanthro$height==low.len & wfhanthro$sex==mat$sex[i]] )
          } else {
            l.val<-wfhanthro$l[wfhanthro$height==low.len & wfhanthro$sex==mat$sex[i]]
            m.val<-wfhanthro$m[wfhanthro$height==low.len & wfhanthro$sex==mat$sex[i]]
            s.val<-wfhanthro$s[wfhanthro$height==low.len & wfhanthro$sex==mat$sex[i]]
          }

          mat$zwfl[i]<-(((mat$weight[i]/m.val)^l.val)-1)/(s.val*l.val)
          if(!is.na(mat$zwfl[i]) & mat$zwfl[i]>3) {
            sd3pos<- m.val*((1+l.val*s.val*3)^(1/l.val))
            sd23pos<- sd3pos- m.val*((1+l.val*s.val*2)^(1/l.val))
            mat$zwfl[i]<- 3+((mat$weight[i]-sd3pos)/sd23pos)
          }
          if(!is.na(mat$zwfl[i]) & mat$zwfl[i]<(-3)) {
            sd3neg<- m.val*((1+l.val*s.val*(-3))**(1/l.val))
            sd23neg<- m.val*((1+l.val*s.val*(-2))**(1/l.val))-sd3neg
            mat$zwfl[i]<- (-3)-((sd3neg-mat$weight[i])/sd23neg)
          }
        }
      }
    }

    if(!is.na(mat$age.day[i]) & mat$age.days[i]>1856) mat$zwfl[i]<-NA

  }

  return(mat)
}

#'  calc.zbmi
#'
#'  Function for calulating individual BMI-for-age z-scores
#' @keywords zscore igrowup who
#' @details WHO Child Growth Standards Department of Nutrition for Health and Development Last modified on 07/10/2013-Developed using R version 3.0.1
#' @note This code conrcerns the standard approach for the prevalences, i.e. the calculation of the prevalences takes into account all the valid (non-missing) z-scores for each of the indicators.
#' @export

calc.zbmi<-function(mat,bmianthro){

  for(i in 1:length(mat$age.days)) {

    if(!is.na(mat$age.days[i]) & mat$age.days[i]>=0 & mat$age.days[i]<=1856 & mat$oedema[i]!="y") {

      l.val<-bmianthro$l[bmianthro$age==mat$age.days[i] & bmianthro$sex==mat$sex[i]]
      m.val<-bmianthro$m[bmianthro$age==mat$age.days[i] & bmianthro$sex==mat$sex[i]]
      s.val<-bmianthro$s[bmianthro$age==mat$age.days[i] & bmianthro$sex==mat$sex[i]]

      mat$zbmi[i]<-(((mat$cbmi[i]/m.val)^l.val)-1)/(s.val*l.val)
      if(!is.na(mat$zbmi[i]) & mat$zbmi[i]>3) {
        sd3pos<- m.val*((1+l.val*s.val*3)^(1/l.val))
        sd23pos<- sd3pos- m.val*((1+l.val*s.val*2)^(1/l.val))
        mat$zbmi[i]<- 3+((mat$cbmi[i]-sd3pos)/sd23pos)
      }
      if(!is.na(mat$zbmi[i]) & mat$zbmi[i]< (-3)) {
        sd3neg<- m.val*((1+l.val*s.val*(-3))**(1/l.val))
        sd23neg<- m.val*((1+l.val*s.val*(-2))**(1/l.val))-sd3neg
        mat$zbmi[i]<- (-3)+((mat$cbmi[i]-sd3neg)/sd23neg)
      }

    } else mat$zbmi[i]<-NA

  }

  return(mat)
}

#'  Modified igrowup
#'
#'  Calculate the z-scores for the indicators: length/height-for-age, weight-for-age, weight-for-legnth/height and body mass index-for-age
#' @keywords zscore igrowup who
#' @details WHO Child Growth Standards Department of Nutrition for Health and Development Last modified on 07/10/2013-Developed using R version 3.0.1
#' @note This code conrcerns the standard approach for the prevalences, i.e. the calculation of the prevalences takes into account all the valid (non-missing) z-scores for each of the indicators.
#' @export

igstd <- function(mydf,sex,age,age.month=F,weight=rep(NA,dim(mydf)[1]),lenhei=rep(NA,dim(mydf)[1]),measure=rep(NA,dim(mydf)[1]),
                             headc=rep(NA,dim(mydf)[1]),armc=rep(NA,dim(mydf)[1]),triskin=rep(NA,dim(mydf)[1]),subskin=rep(NA,dim(mydf)[1]),oedema=rep("n",dim(mydf)[1]),sw=rep(1,dim(mydf)[1])) {
#Calculating the z-scores for all indicators
  old <- options(warn=(-1))

  sex.x<-as.character(get(deparse(substitute(mydf)))[,deparse(substitute(sex))])
  age.x<-as.double(get(deparse(substitute(mydf)))[,deparse(substitute(age))])
  if(!missing(weight)) weight.x<-as.double(get(deparse(substitute(mydf)))[,deparse(substitute(weight))]) else weight.x<-as.double(weight)
  if(!missing(lenhei)) lenhei.x<-as.double(get(deparse(substitute(mydf)))[,deparse(substitute(lenhei))]) else lenhei.x<-as.double(lenhei)
  if(!missing(headc)) headc.x<-as.double(get(deparse(substitute(mydf)))[,deparse(substitute(headc))]) else headc.x<-as.double(headc)
  if(!missing(armc)) armc.x<-as.double(get(deparse(substitute(mydf)))[,deparse(substitute(armc))]) else armc.x<-as.double(armc)
  if(!missing(triskin)) triskin.x<-as.double(get(deparse(substitute(mydf)))[,deparse(substitute(triskin))]) else triskin.x<-as.double(triskin)
  if(!missing(subskin)) subskin.x<-as.double(get(deparse(substitute(mydf)))[,deparse(substitute(subskin))]) else subskin.x<-as.double(subskin)
  if(!missing(measure)) lorh.vec<-as.character(get(deparse(substitute(mydf)))[,deparse(substitute(measure))]) else lorh.vec<-as.character(measure)
  if(!missing(oedema)) oedema.vec<-as.character(get(deparse(substitute(mydf)))[,deparse(substitute(oedema))]) else oedema.vec<-oedema
  if(!missing(sw))	sw<-as.double(get(deparse(substitute(mydf)))[,deparse(substitute(sw))])	else sw<-as.double(sw)
  sw<-ifelse(is.na(sw),0,sw)

  sex.vec<-NULL

  if(age.month) age.vec<-rounde(age.x*30.4375) else age.vec<-rounde(age.x)
  lenhei.vec<-ifelse((!is.na(age.vec) & age.vec<731 & !is.na(lorh.vec) & (lorh.vec=="h" | lorh.vec=="H")),lenhei.x+0.7,#
                     ifelse((!is.na(age.vec) & age.vec>=731 & !is.na(lorh.vec) & (lorh.vec=="l" | lorh.vec=="L")),lenhei.x-0.7,lenhei.x))

  sex.vec<-ifelse(!is.na(sex.x) & (sex.x=="m" | sex.x=="M" | sex.x=="1"),1,ifelse(!is.na(sex.x) & (sex.x=="f" | sex.x=="F" | sex.x=="2"),2,NA))

  lorh.vec<-ifelse(is.na(lorh.vec) | lorh.vec=="l" | lorh.vec=="L" | lorh.vec=="h" | lorh.vec=="H",lorh.vec,NA)

  oedema.vec<-ifelse(oedema.vec=="n" | oedema.vec=="N","n",ifelse(oedema.vec=="y" | oedema.vec=="Y","y","n"))

  mat<-cbind.data.frame(age.x,as.integer(age.vec),as.double(sex.vec),weight.x,lenhei.x,lorh.vec,lenhei.vec,headc.x,armc.x,triskin.x,subskin.x,oedema.vec,sw,stringsAsFactors=F)
  names(mat)<-c("age","age.days","sex","weight","len.hei","l.h","clenhei","headc","armc","triskin","subskin","oedema","sw")

  mat$cbmi<-mat$weight/((lenhei.vec/100)^2)
  mat$zlen<-rep(NA,length(mat$age))
  mat$zwei<-rep(NA,length(mat$age))
  mat$zwfl<-rep(NA,length(mat$age))
  mat$zbmi<-rep(NA,length(mat$age))
  mat$zhc<-rep(NA,length(mat$age))
  mat$zac<-rep(NA,length(mat$age))
  mat$zts<-rep(NA,length(mat$age))
  mat$zss<-rep(NA,length(mat$age))

  mat$zlen<-rep(NA,length(mat$age))
  mat$flen<-rep(NA,length(mat$age))
  mat$fwei<-rep(NA,length(mat$age))
  mat$fwfl<-rep(NA,length(mat$age))
  mat$fbmi<-rep(NA,length(mat$age))
  mat$fhc<-rep(NA,length(mat$age))
  mat$fac<-rep(NA,length(mat$age))
  mat$fts<-rep(NA,length(mat$age))
  mat$fss<-rep(NA,length(mat$age))

#Calculating the z-scores for all indicators
  cat("Please wait while calculating z-scores...\n")
  ### Length-for-age z-score
  mat<-calc.zlen(mat,lenanthro)
  ### Head circumference-for-age z-score
  mat<-calc.zhc(mat,hcanthro)
  ### Weight-for-age z-score
  mat<-calc.zwei(mat,weianthro)
  ### Arm circumference-for-age z-score
  mat<-calc.zac(mat,acanthro)
  ### Triceps skinfold-for-age z-score
  mat<-calc.zts(mat,tsanthro)
  ### Subscapular skinfold-for-age z-score
  mat<-calc.zss(mat,ssanthro)
  ### Weight-for-length/height z-score
  mat<-calc.zwfl(mat,wflanthro,wfhanthro)
  ### BMI-for-age z-score
  mat<-calc.zbmi(mat,bmianthro)
  #### Rounding the z-scores to two decimals
  mat$zlen<-rounde(mat$zlen,digits=2)
  mat$zwei<-rounde(mat$zwei,digits=2)
  mat$zwfl<-rounde(mat$zwfl,digits=2)
  mat$zbmi<-rounde(mat$zbmi,digits=2)
  mat$zhc<-rounde(mat$zhc,digits=2)
  mat$zac<-rounde(mat$zac,digits=2)
  mat$zts<-rounde(mat$zts,digits=2)
  mat$zss<-rounde(mat$zss,digits=2)

  #### Flagging z-score values for individual indicators
  mat$flen<-ifelse(abs(mat$zlen) > 6,1,0)
  mat$fwei<-ifelse(mat$zwei > 5 | mat$zwei < (-6),1,0)
  mat$fwfl<-ifelse(abs(mat$zwfl) > 5,1,0)
  mat$fbmi<-ifelse(abs(mat$zbmi) > 5,1,0)
  mat$fhc<-ifelse(abs(mat$zhc) > 5,1,0)
  mat$fac<-ifelse(abs(mat$zac) > 5,1,0)
  mat$fts<-ifelse(abs(mat$zts) > 5,1,0)
  mat$fss<-ifelse(abs(mat$zss) > 5,1,0)

  cat(paste("Z-scores calculated - preparing matrix"))
  mat<-cbind.data.frame(mydf,mat[,-c(1,3:6,8:9)])

  ######### Export data frame with z-scores and flag variables
  assign("matz",mat,envir = .GlobalEnv)
}   #### End of function igrowup.standard

#' Standard data igrowup
#'
#' Data from the standard igrowup files
#'
#' @docType data
#'
#' @usage data(igdata)
#'
#' @format An object of class \code{"cross"}; see \code{\link[qtl]{read.cross}}.
#'
#'
#' @keywords datasets
#'
#' @source \href{https://www.who.int/childgrowth/software/igrowup_R.zip}{who}
#'
#' @examples
#' data("igdata")

