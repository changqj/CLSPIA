function para=Parameterize(data,class)%Parameterize given npoints (x,y). class: 0:uniform  0.5:centr    1.0:chord
l_sum=0; %sum of all the chord
n=numel(data(:,1));
para([1 n])=0; %parameters of the points
l([1 n])=0;
for i=2:1:n
    l(i-1)=norm(data(i,:)-data(i-1,:))^(class);
    l_sum=l_sum+l(i-1); %compute the sum of all chords
end
for i=2:n-1
    para(i)=para(i-1)+l(i-1)/l_sum; %compute parameters for all points     
end
para(n) = 1;