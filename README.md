# R-code
#winsorized mean

w.mean=function(x,g)#g is winsorized location in %
{
n=length(x)
y=sort(x)
wsm=((g+1)*y[(g+1)]+sum(y[((g+2):(n-g-1))])+(g+1)*y[(n-g)])/n
wsm
}

#winsorized variance

w.ss=function(x,g){
y=sort(x)
n=length(x)
z=w.mean(x,g)
k=function(x,g){t=0;
for(i in (g+2):(n-g-1)) t=t+(y[(i)]-z)^2
t
}
wd=((g+1)*(y[(g+1)]-z)^2+k(x,g)+(g+1)*(y[(n-g)]-z)^2)/(n-2*g-1)
wd
}


#population variance uw
uw=function(x,g){
y=sort(x)
n=length(x)
z=w.mean(x,g)
k=function(x,g){t=0;
for(i in (g+2):(n-g-1)) t=t+(y[(i)]-z)^3
t
}
wd=((g+1)*(y[(g+1)]-z)^3+k(x,g)+(g+1)*(y[(n-g)]-z)^3)/(n-2*g)
wd
}


#degree of freedom v
fd=function(x,y,g1,g2){
n1=length(x)
n2=length(y)
x1=mean(x,trim=g1/n1)
y1=mean(y,trim=g2/n2)
f1=n1-2*g1
f2=n2-2*g2
s1=w.ss(x,g1)/f1
s2=w.ss(y,g2)/f2
ns=s1+s2
v1=(s1)^2
v2=(s2)^2
v=(ns*ns)/{(v1/(f1-1))+(v2/(f2-1))}
v

}



#Johnson transformation t test
johnson.t=function(x,y,g1,g2){
n1=length(x)
n2=length(y)
x1=mean(x,trim=g1/n1)
y1=mean(y,trim=g2/n2)
f1=n1-2*g1
f2=n2-2*g2
s1=w.ss(x,g1)/f1
s2=w.ss(y,g2)/f2
uw1=uw(x,g1)
uw2=uw(y,g2)
nuw=(uw1/(f1)^2)-(uw2/(f2)^2)
ns=s1+s2
k=nuw/(6*ns)
s=nuw/(3*ns*ns)
ht=((x1-y1)+k+((x1-y1)^2*s))/sqrt(ns)
ht
}

#Johnson transformation t test with pool standard deviation
johnson.t.1=function(x,y,g1,g2){
n1=length(x)
n2=length(y)
x1=mean(x,trim=g1/n1)
y1=mean(y,trim=g2/n2)
f1=n1-2*g1
f2=n2-2*g2
s1=w.ss(x,g1)/f1
s2=w.ss(y,g2)/f2
sp2=(w.ss(x,g1)*(n1-1)+w.ss(y,g1)*(n2-1))/(n1+n2-2)
uw1=uw(x,g1)
uw2=uw(y,g2)
nuw=(uw1/(f1)^2)-(uw2/(f2)^2)
ns=sp2*(1/f1+1/f2)
k=nuw/(6*ns)
s=nuw/(3*ns*ns)
ht=((x1-y1)+k+((x1-y1)^2*s))/sqrt(ns)
ht
}


#Johnson transformation t test with control standard deviation
johnson.t.2=function(x,y,g1,g2){
n1=length(x)
n2=length(y)
x1=mean(x,trim=g1/n1)
y1=mean(y,trim=g2/n2)
f1=n1-2*g1
f2=n2-2*g2
s1=w.ss(x,g1)/f1
s2=w.ss(y,g2)/f2
sp2=w.ss(x,g1) #or sp2=w.ss(y,g2)
uw1=uw(x,g1)
uw2=uw(y,g2)
nuw=(uw1/(f1)^2)-(uw2/(f2)^2)
ns=sp2*(1/f1+1/f2)
k=nuw/(6*ns)
s=nuw/(3*ns*ns)
ht=((x1-y1)+k+((x1-y1)^2*s))/sqrt(ns)
ht
}


#half student t test 
half.t=function(x,y,g1,g2){
n1=length(x)
n2=length(y)
x1=mean(x,trim=g1/n1)
y1=mean(y,trim=g2/n2)
f1=n1-2*g1
f2=n2-2*g2
s1=w.ss(x,g1)/f1
s2=w.ss(y,g2)/f2
sp2=w.ss(x,g1) #or sp2=w.ss(y,g2)
uw1=uw(x,g1)
uw2=uw(y,g2)
nuw=(uw1/(f1)^2)-(uw2/(f2)^2)
ns=sp2*(1/f1+1/f2)
k=nuw/(6*ns)
s=nuw/(3*ns*ns)
ht=(x1-y1)/sqrt(ns)
ht
}


#pool student t test 
pool.t=function(x,y,g1,g2){
n1=length(x)
n2=length(y)
x1=mean(x,trim=g1/n1)
y1=mean(y,trim=g2/n2)
f1=n1-2*g1
f2=n2-2*g2
s1=w.ss(x,g1)/f1
s2=w.ss(y,g2)/f2
sp2=(w.ss(x,g1)*(n1-1)+w.ss(y,g1)*(n2-1))/(n1+n2-2)
uw1=uw(x,g1)
uw2=uw(y,g2)
nuw=(uw1/(f1)^2)-(uw2/(f2)^2)
ns=sp2*(1/f1+1/f2)
k=nuw/(6*ns)
s=nuw/(3*ns*ns)
ht=(x1-y1)/sqrt(ns)
ht
}

#normal for type 1 errpr and power
tmpa=0
tmpa1=0
tmpa2=0
tmpa3=0
tmpb=0
tmpb1=0
tmpb2=0
tmpb3=0
num1=20
num2=20
r=1 
for(i in 1:100000)
{

y0=rnorm(num1,100,30)
y1=rnorm(num2,100,30*r)


a=(1-pt(abs(johnson.t.2(y0,y1,0,0)), num1-1)<0.05/2)^1

tmpa=tmpa+a


a1=(1-pt(abs(johnson.t.2(y0,y1,0,0)), num1-1)<0.01/2)^1

tmpa1=tmpa1+a1


a2=(1-pt(abs(johnson.t.2(y0,y1,0,0)), num1-1)<0.005/2)^1

tmpa2=tmpa2+a2


a3=(1-pt(abs(johnson.t.2(y0,y1,0,0)), num1-1)<0.001/2)^1

tmpa3=tmpa3+a3


b=(1-pt(abs(half.t(y0,y1,0,0)), num1-1)<0.05/2)^1
tmpb=tmpb+b


b1=(1-pt(abs(half.t(y0,y1,0,0)), num1-1)<0.01/2)^1
tmpb1=tmpb1+b1



b2=(1-pt(abs(half.t(y0,y1,0,0)), num1-1)<0.005/2)^1
tmpb2=tmpb2+b2


b3=(1-pt(abs(half.t(y0,y1,0,0)), num1-1)<0.001/2)^1
tmpb3=tmpb3+b3

}

sink("d:/type1_n_(1).txt")
c(tmpa,tmpa1,tmpa2,tmpa3,tmpb,tmpb1,tmpb2,tmpb3)*100/100000
sink()

#uniform for type 1 errpr and power
tmpa=0
tmpa1=0
tmpa2=0
tmpa3=0
tmpb=0
tmpb1=0
tmpb2=0
tmpb3=0
num1=60
num2=20
a0=48.03847577
b0=103.9230485+48.03847577
a1=48.03847577
b1=103.9230485+48.03847577
for(i in 1:100000)
{

y0=runif(num1,a0,b0)
y1=runif(num2,a1,b1)


a=(1-pt(abs(johnson.t.2(y0,y1,0,0)), num1-1)<0.05/2)^1

tmpa=tmpa+a


a1=(1-pt(abs(johnson.t.2(y0,y1,0,0)), num1-1)<0.01/2)^1

tmpa1=tmpa1+a1


a2=(1-pt(abs(johnson.t.2(y0,y1,0,0)), num1-1)<0.005/2)^1

tmpa2=tmpa2+a2


a3=(1-pt(abs(johnson.t.2(y0,y1,0,0)), num1-1)<0.001/2)^1

tmpa3=tmpa3+a3


b=(1-pt(abs(half.t(y0,y1,0,0)), num1-1)<0.05/2)^1
tmpb=tmpb+b


b1=(1-pt(abs(half.t(y0,y1,0,0)), num1-1)<0.01/2)^1
tmpb1=tmpb1+b1



b2=(1-pt(abs(half.t(y0,y1,0,0)), num1-1)<0.005/2)^1
tmpb2=tmpb2+b2


b3=(1-pt(abs(half.t(y0,y1,0,0)), num1-1)<0.001/2)^1
tmpb3=tmpb3+b3

}


sink("d:/type_1_u_(1).txt")
c(tmpa,tmpa1,tmpa2,tmpa3,tmpb,tmpb1,tmpb2,tmpb3)*100/100000
sink()

#gamma for type 1 error and power
tmpa=0
tmpa1=0
tmpa2=0
tmpa3=0
tmpb=0
tmpb1=0
tmpb2=0
tmpb3=0

num1=60
num2=20
for(i in 1:100000)
{

y0=rgamma(num1,(100^2)/(30*1)^2,1/((30*1)^2/100))
y1=rgamma(num2,(100^2)/(30*1)^2,1/((30*1)^2/100))


a=(1-pt(abs(johnson.t.2(y0,y1,0,0)), num1-1)<0.05/2)^1

tmpa=tmpa+a


a1=(1-pt(abs(johnson.t.2(y0,y1,0,0)), num1-1)<0.01/2)^1

tmpa1=tmpa1+a1


a2=(1-pt(abs(johnson.t.2(y0,y1,0,0)), num1-1)<0.005/2)^1

tmpa2=tmpa2+a2


a3=(1-pt(abs(johnson.t.2(y0,y1,0,0)), num1-1)<0.001/2)^1

tmpa3=tmpa3+a3


b=(1-pt(abs(half.t(y0,y1,0,0)), num1-1)<0.05/2)^1
tmpb=tmpb+b


b1=(1-pt(abs(half.t(y0,y1,0,0)), num1-1)<0.01/2)^1
tmpb1=tmpb1+b1



b2=(1-pt(abs(half.t(y0,y1,0,0)), num1-1)<0.005/2)^1
tmpb2=tmpb2+b2


b3=(1-pt(abs(half.t(y0,y1,0,0)), num1-1)<0.001/2)^1
tmpb3=tmpb3+b3



}


sink("d:/type1_g_(1).txt")
c(tmpa,tmpa1,tmpa2,tmpa3,tmpb,tmpb1,tmpb2,tmpb3)*100/100000
sink()

#negative gamma for type 1 error and power
tmpa=0
tmpa1=0
tmpa2=0
tmpa3=0
tmpb=0
tmpb1=0
tmpb2=0
tmpb3=0
num1=20
num2=20
r=1 
for(i in 1:100000)
{

y0=200-rgamma(num1,(100^2)/(30*1)^2,1/((30*1)^2/100))
y1=200-rgamma(num2,(100^2)/(30*1*r)^2,1/((30*1*r)^2/100))



a=(1-pt(abs(johnson.t.2(y0,y1,0,0)), num1-1)<0.05/2)^1

tmpa=tmpa+a


a1=(1-pt(abs(johnson.t.2(y0,y1,0,0)), num1-1)<0.01/2)^1

tmpa1=tmpa1+a1


a2=(1-pt(abs(johnson.t.2(y0,y1,0,0)), num1-1)<0.005/2)^1

tmpa2=tmpa2+a2


a3=(1-pt(abs(johnson.t.2(y0,y1,0,0)), num1-1)<0.001/2)^1

tmpa3=tmpa3+a3


b=(1-pt(abs(half.t(y0,y1,0,0)), num1-1)<0.05/2)^1
tmpb=tmpb+b


b1=(1-pt(abs(half.t(y0,y1,0,0)), num1-1)<0.01/2)^1
tmpb1=tmpb1+b1



b2=(1-pt(abs(half.t(y0,y1,0,0)), num1-1)<0.005/2)^1
tmpb2=tmpb2+b2


b3=(1-pt(abs(half.t(y0,y1,0,0)), num1-1)<0.001/2)^1
tmpb3=tmpb3+b3

}


sink("f:/type1_ng_(1).txt")
c(tmpa,tmpa1,tmpa2,tmpa3,tmpb,tmpb1,tmpb2,tmpb3)*100/100000
sink()
