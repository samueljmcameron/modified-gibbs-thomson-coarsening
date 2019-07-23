from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib


xl = -2
xh = 12
yl = -0.1
yh = 0.1
nm_points = 20
#ss = 20. # 1/8
ss = 0.01
Req = 6.
d = 2.0

R = np.linspace(0.1,xh,num = 300)
"""
def R_u(Req,ss,d):
    a = (d+1.0)/(d+2.0)
    return Req*(a-np.sqrt(a**2-d/(d+2.0)*(1-2*ss/(Req**2*gppReq))))

def R_s(Req,ss,d):
    a = (d+1.0)/(d+2.0)
    return Req*(a+np.sqrt(a**2-d/(d+2.0)*(1-2*ss/(Req**2*gppReq))))
"""
def Rdot(R,Req,ss,d,alpha):
    return 1.0/R*(ss-(1-1.0/d)/R
                  +(1.0/alpha-1.0/d)/(Req**(1-alpha)*R**alpha)
                  +(1-1.0/alpha)/Req)

#plt.xkcd()

fig = plt.figure()
ax = fig.add_subplot(1,1,1)

ax.spines['right'].set_color('none')
ax.spines['left'].set_color('none')
ax.spines['top'].set_color('none')
ax.spines['bottom'].set_color('none')
ax.set_xticks([])
ax.set_yticks([])
ax.set_ylim([yl,yh])
ax.set_xlim([xl,xh])


xaxis = np.linspace(xl,xh,num=nm_points,endpoint = True)
yaxis = np.linspace(yl,yh,num=nm_points,endpoint = True)
lining = np.ones(nm_points)

ax.plot(xaxis,lining*0,'k',lw=2.0)
ax.plot(lining*0,yaxis,'k',lw=2.0)
ax.text(xl+0.01,yh-0.01,r'$\dot{R}$',fontsize=25)
ax.text(xh-1,-0.02,r'$R$',fontsize=25)

"""
ax.text(R_u(Req,ss,d)-1,-.35,'<',fontsize=30)
ax.text(R_u(Req,ss,d)-2,-.35,'<',fontsize=30)
ax.text(R_u(Req,ss,d)+0.75,-.35,'>',fontsize=30)
ax.text(R_u(Req,ss,d)+1.75,-.35,'>',fontsize=30)
ax.annotate('unstable',xy=(R_u(Req,ss,d)-0.1,+0.1),xytext=(R_u(Req,ss,d)-2.2,+2),
            arrowprops=dict(facecolor='black',shrink=0.05))

ax.text(R_s(Req,ss,d)-1,-.35,'>',fontsize=30)
ax.text(R_s(Req,ss,d)-2,-.35,'>',fontsize=30)
ax.text(R_s(Req,ss,d)+0.75,-.35,'<',fontsize=30)
ax.text(R_s(Req,ss,d)+1.75,-.35,'<',fontsize=30)
ax.annotate('stable',xy=(R_s(Req,ss,d)+0.1,+0.1),xytext=(R_s(Req,ss,d)+1,+3),
            arrowprops=dict(facecolor='black',shrink=0.05))
"""

ss1 = 0.01
ax.plot(R,Rdot(R,Req,ss1,d,1.0/3.0),'b',lw=3.0,
        label=r"$\chi(t)=%.2lf$"%ss1)
ss2 = 0.1
ax.plot(R,Rdot(R,Req,ss2,d,1.0/3.0),'m',lw=3.0,
        label=r"$\chi(t)=%.2lf$"%ss2)
ax.legend(frameon=False,fontsize=18)
#ax.plot(R_u(Req,ss,d),0,'ko')
#ax.text(R_u(Req,ss,d),-1,r'$R_u$',fontsize=20)
#ax.plot(R_s(Req,ss,d),0,'ko')
#ax.text(R_s(Req,ss,d)-0.4,-1,r'$R_s$',fontsize=20)
#ax.plot(Req,0,'ko')#,markersize=20,markeredgewidth=2.)
#ax.text(Req,-1,r'$R_{eq}$',fontsize = 20)


fig.savefig('fitted_collagen.pdf')
plt.show()

