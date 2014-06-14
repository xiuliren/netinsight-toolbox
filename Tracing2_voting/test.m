clc

subplot(211)
a=gca() %get the current axes
a.box='off';
t=-%pi:0.3:%pi;plot3d(t,t,sin(t)'*cos(t),80,50,'X@Y@Z',[5,2,4]);
subplot(212)
plot2d(); %simple plot
a=gca() %get the current axes
a.box='off';
a.x_location='middle';
a.parent.background=4;
delete(gca()) %%// delete the current axes    
xdel(0) %delete a graphics window