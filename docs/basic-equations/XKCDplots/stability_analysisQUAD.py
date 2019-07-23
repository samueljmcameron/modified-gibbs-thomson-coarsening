from __future__ import print_function, division
import numpy as np
import matplotlib.pyplot as plt
import matplotlib


xl = -2
xh = 12 
yl = -7
yh = 7
nm_points = xh-xl+1
#ss = 20. # 1/8
ss = 0.
Req = 6.
d = 3.0
gppReq = 4. #20.

R = np.linspace(0.1,xh,num = 300)

def R_u(Req,ss,d):
    a = (d+1.0)/(d+2.0)
    return Req*(a-np.sqrt(a**2-d/(d+2.0)*(1-2*ss/(Req**2*gppReq))))

def R_s(Req,ss,d):
    a = (d+1.0)/(d+2.0)
    return Req*(a+np.sqrt(a**2-d/(d+2.0)*(1-2*ss/(Req**2*gppReq))))

def Rdot(R,Req,ss,d):
    return 1.0/R*(ss-gppReq/2.0*((d+2.0)/d*R-Req)*(R-Req))

print(R_u(Req,ss,d))

plt.xkcd()

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

ax.plot(xaxis,lining*0,'k')
ax.plot(lining*0,yaxis,'k')
ax.text(xl+1,yh-1,r'$\dot{R}$',fontsize=25)
ax.text(xh-1,yl+6,r'$R$',fontsize=25)

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

ax.plot(R,Rdot(R,Req,ss,d))
ax.plot(R_u(Req,ss,d),0,'ko')
ax.text(R_u(Req,ss,d),-1,r'$R_u$',fontsize=20)
#ax.plot(R_s(Req,ss,d),0,'ko')
#ax.text(R_s(Req,ss,d)-0.4,-1,r'$R_s$',fontsize=20)
ax.plot(Req,0,'ko')#,markersize=20,markeredgewidth=2.)
ax.text(Req,-1,r'$R_{eq}$',fontsize = 20)


fig.savefig('stability_at_eq.pdf')
plt.show()

