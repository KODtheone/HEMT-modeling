%2016年1月16日16:17:28  完成了，根据 外部寄生参数 和 hot-fet测试S参数， 计算内部寄生参数的值
%2016年3月3日20:57:10 修改程序，获得外部寄生参数的 优化提取值 命名：main_V2
%2016年3月9日22:29:44 verify版本 用于验证使用M函数后的程序正确性   验证！
%2016年3月17日22:43:02 用于加入扫描 RTH
%2016年4月21日14:53:05  采用新的top  尝试调节权重以获得合适的 变化温度?T
%2016年4月23日20:17:22  版本：progress_mainV3_intYS   表示在int参数取值时，误差考虑的是Y和S加起来的。
%2016年4月28日19:47:57 统一ε1和ε2的计算方法；
%2017年5月17日21:40:45 内部参数直接使用部分平均值

S_test = sparameters('D:\学习\MATLAB\2017年4月18日 论文SP结果 最简单程序\VGS = -1.5VDS = 5.txt');
T=25;
Rdd=6.03+0.0614*T
Rss=1.99+0.0142*T    %扫描的是T  #2016年5月13日11:13:52，注意基准温度是0
% Rdd=7.56;
% Rss=2.34;
% S_test = sparameters('VGSn2VDS10.s2p');
%完成 test数据导入
Rd = Rdd;
Ld= 58.1e-12;
Cpd= 83e-15;
%testABCD=s2abcd(SD,50);
%A = magic(5)
%F=SD.Parameters
%testABCD=s2abcd(F,50)
%K=F.*testABCD.*F
Rg =1.6;
Lg=46.9e-12;
Cpg=50e-15; ;
%G电路
Rs = Rss;
Ls=2.79e-12;
%S电路
%2017年5月17日19:57:10， 改为直接Z参数计算
freq=linspace(0.1e9,20.1e9,51);
w=2*pi*freq;
Z22d=j.*w.*Ld+1./(j.*w.*Cpd);
Z12d=1./(j.*w.*Cpd);
Z21d=1./(j.*w.*Cpd);
Z11d=Rd+1./(j.*w.*Cpd);
Z11dd=Z11d.';
Z12dd=Z12d.';
Z21dd=Z21d.';
Z22dd=Z22d.';
Dz=zeros(2,2,51);
for k=1:51
  Dz(1,1,k) =Z11dd(k,1);
  Dz(1,2,k) =Z12dd(k,1);
  Dz(2,1,k) =Z21dd(k,1);
  Dz(2,2,k) =Z22dd(k,1);
end
Dabcd=z2abcd(Dz);
% D 电路

Z11g=j.*w.*Lg+1./(j.*w.*Cpg);
Z12g=1./(j.*w.*Cpg);
Z21g=1./(j.*w.*Cpg);
Z22g=Rg+1./(j.*w.*Cpg);
Z11gg=Z11g.';
Z12gg=Z12g.';
Z21gg=Z21g.';
Z22gg=Z22g.';
Gz=zeros(2,2,51);
for k=1:51
  Gz(1,1,k) =Z11gg(k,1);
  Gz(1,2,k) =Z12gg(k,1);
  Gz(2,1,k) =Z21gg(k,1);
  Gz(2,2,k) =Z22gg(k,1);
end
Gabcd=z2abcd(Gz);
%G电路
Z11s=Rs+j.*w.*Ls;
Z12s=Rs+j.*w.*Ls;
Z21s=Rs+j.*w.*Ls;
Z22s=Rs+j.*w.*Ls;

Z11ss=Z11s.';
Z12ss=Z12s.';
Z21ss=Z21s.';
Z22ss=Z22s.';
Sz=zeros(2,2,51);
for k=1:51
  Sz(1,1,k) =Z11ss(k,1);
  Sz(1,2,k) =Z12ss(k,1);
  Sz(2,1,k) =Z21ss(k,1);
  Sz(2,2,k) =Z22ss(k,1);
end
% S电路
ALLabcd=s2abcd(S_test.Parameters);

IN_D(2,2,51)=0;
IN_G(2,2,51)=0;
ISabcd(2,2,51)=0;
for k=1:51
  IN_D(:,:,k)=[1,0;0,1]/Dabcd(:,:,k);
  IN_G(:,:,k)=[1,0;0,1]/Gabcd(:,:,k);
  ISabcd(:,:,k)=IN_G(:,:,k)*ALLabcd(:,:,k)*IN_D(:,:,k);
end

ISz=abcd2z(ISabcd);
IFETz=ISz-Sz;
IFETy=z2y(IFETz);   %矩阵计算

re011=real(IFETy(1,1,:));
re11=squeeze(re011(1,:,:))';
re012=real(IFETy(1,2,:));
re12=squeeze(re012(1,:,:))';
re021=real(IFETy(2,1,:));
re21=squeeze(re021(1,:,:))';
re022=real(IFETy(2,2,:));
re22=squeeze(re022(1,:,:))';
im011=imag(IFETy(1,1,:));
im11=squeeze(im011(1,:,:))';
im012=imag(IFETy(1,2,:));
im12=squeeze(im012(1,:,:))';
im021=imag(IFETy(2,1,:));
im21=squeeze(im021(1,:,:))';
im022=imag(IFETy(2,2,:));
im22=squeeze(im022(1,:,:))';
% ===参数计算部分；
Cgd=(-im12)./w.*(1+re12.*re12./im12./im12);
Cgs=(power(re11+re12,2)+power(im11+im12,2))./w./(im11+im12);
Rds=1./(re22+re12);
Cds=(im22+im12)./w;
Ri=(re11+re12)./(power(im11+im12,2)+power(re11+re12,2));
gm=sqrt((power(re21-re12,2)-power(im21-im12,2)).*(1+power(w.*Cgs.*Ri,2)));
t=-1./w.*atan((im21-im12+w.*Ri.*Cgs.*(re21-re12))./(re21-re12-w.*Ri.*Cgs.*(im21-im12)));
Rgd=-re12./(power(im12,2).*(1+power(re12./im12,2)));
%参数提取完成
Cgs_mean=mean(Cgs(12:42));
Cds_mean=mean(Cds(12:42));
t_mean=mean(t(12:42));
Cgd_mean=mean(Cgd(12:42));
gm_mean=mean(gm(12:42));
Ri_mean=mean(Ri(12:42));
Rds_mean=mean(Rds(12:42));
Rgd_mean=mean(Rgd(12:42));

%制造函数：



RE=100;
yRE=100;
Re=100;
allabcd(2,2,51)=0;
% h = waitbar(0,'Please wait...');
%     
for i1=1:n(1)
    for i2=1:n(2)
       
        for i3=1:n(3)
            
           for i4=1:n(4)
                
               for i5=1:n(5)
%                waitbar(((i1-1)*n(2)*n(3)*n(4)*n(5)+(i2-1)*n(3)*n(4)*n(5)+(i3-1)*n(4)*n(5)+(i4-1)*n(5))/(n(1)*n(2)*n(3)*n(4)*n(5)),h);
                     for i6=1:n(6)
                         for i7=1:n(7)
                             for i8=1:n(8)
                             var=[Cgs_n(i1),Cds_n(i2),t_n(i3),Cgd_n(i4),gm_n(i5),Ri_n(i6),Rds_n(i7),Rgd_n(i8)];
                              %RE=calEs(var);   2016年3月19日14:48:27  将这句话用程序替换；
                      
%                        ==============  intrinsic 反向计算
                              Y110=1i.*w.*Cgs_n(i1)./(1+1i.*w.*Cgs_n(i1).*Ri_n(i6))+1i.*w.*Cgd_n(i4)./(1+1i.*w.*Cgd_n(i4).*Rgd_n(i8));
                              Y120=-1i.*w.*Cgd_n(i4)./(1+1i.*w.*Cgd_n(i4).*Rgd_n(i8));
                              Y210=gm_n(i5).*exp(-1i.*w.*t_n(i3))./(1+1i.*w.*Cgs_n(i1).*Ri_n(i6))-1i.*w.*Cgd_n(i4)./(1+1i.*w.*Cgd_n(i4).*Rgd_n(i8));
                              Y220=1./Rds_n(i7)+1i.*w.*Cds_n(i2)+1i.*w.*Cgd_n(i4)./(1+1i.*w.*Cgd_n(i4).*Rgd_n(i8));
                        %=========intrinsic反向计算结束
                              Y11=Y110.';
                              Y12=Y120.';
                              Y21=Y210.';
                              Y2122=Y21.';
                              Y22=Y220.';
                              ifety=zeros(2,2,51);
                                for k=1:51
                                ifety(1,1,k) =Y11(k,1);
                                ifety(1,2,k) =Y12(k,1);
                                ifety(2,1,k) =Y21(k,1);
                                ifety(2,2,k) =Y22(k,1);
                                end
                              ifetz=y2z(ifety);
                              isz=ifetz+SSz;
                              isabcd=z2abcd(isz);

                                for k=1:51
                                  allabcd(:,:,k)=Gabcd(:,:,k)*isabcd(:,:,k)*Dabcd(:,:,k);
                                end
                              alls=abcd2s(allabcd);
% compare  alls &S_test.Parameters
%============2016年3月23日19:43:16 尝试不归一   SP误差计算  ε2
%误差形式1 
% dS11=alls(1,1,:)-S_test.Parameters(1,1,:);
% ads11=abs(dS11)./abs(S_test.Parameters(1,1,:));
% % ads11=abs(dS11)/max(abs(S_test.Parameters(1,1,:)));
% es11=sum(ads11)/51;
% dS12=alls(1,2,:)-S_test.Parameters(1,2,:);
% ads12=abs(dS12)./abs(S_test.Parameters(1,2,:));
% % ads12=abs(dS12)/max(abs(S_test.Parameters(1,2,:)));
% es12=sum(ads12)/51;
% dS21=alls(2,1,:)-S_test.Parameters(2,1,:);
% ads21=abs(dS21)./abs(S_test.Parameters(2,1,:));
% % ads21=abs(dS21)/max(abs(S_test.Parameters(2,1,:)));
% es21=sum(ads21)/51;
% dS22=alls(2,2,:)-S_test.Parameters(2,2,:);
% ads22=abs(dS22)./abs(S_test.Parameters(2,2,:));
% % ads22=abs(dS22)/max(abs(S_test.Parameters(2,2,:)));
% es22=sum(ads22)/51;
%误差形式2 用 max
                              dS11=alls(1,1,:)-S_test.Parameters(1,1,:);
                              ads11=abs(dS11)/max(abs(S_test.Parameters(1,1,:)));
                              es11=sum(ads11)/51;
                              dS12=alls(1,2,:)-S_test.Parameters(1,2,:);
                              ads12=abs(dS12)/max(abs(S_test.Parameters(1,2,:)));
                              es12=sum(ads12)/51;
                              dS21=alls(2,1,:)-S_test.Parameters(2,1,:);
                              ads21=abs(dS21)/max(abs(S_test.Parameters(2,1,:)));
                              es21=sum(ads21)/51;
                              dS22=alls(2,2,:)-S_test.Parameters(2,2,:);
                              ads22=abs(dS22)/max(abs(S_test.Parameters(2,2,:)));
                              es22=sum(ads22)/51;
%===
                              ws11=0.25;ws12=0.25;ws21=0.25;ws22=0.25;    %2017年4月18日15:57:02 修改S22权重！
                              E=(ws11*es11+ws12*es12+ws21*es21+ws22*es22)*100;%
%==============================2016年3月23日16:56:06 以IFETy的误差为计算标准：     ε1
%误差形式1 
% dFETy11=ifety(1,1,:)-IFETy(1,1,:);
% adFETy11=abs(dFETy11)./abs(IFETy(1,1,:));
% % adFETy11=abs(dFETy11)/max(abs(IFETy(1,1,:)));
% eFETy11=sum(adFETy11)/51;
% dFETy12=ifety(1,2,:)-IFETy(1,2,:);
% adFETy12=abs(dFETy12)./abs(IFETy(1,2,:));
% % adFETy12=abs(dFETy12)/max(abs(IFETy(1,2,:)));
% eFETy12=sum(adFETy12)/51;
% dFETy21=ifety(2,1,:)-IFETy(2,1,:);
% adFETy21=abs(dFETy21)./abs(IFETy(2,1,:));
% % adFETy21=abs(dFETy21)/max(abs(IFETy(2,1,:)));
% eFETy21=sum(adFETy21)/51;
% dFETy22=ifety(2,2,:)-IFETy(2,2,:);
% adFETy22=abs(dFETy22)./abs(IFETy(2,2,:));
% % adFETy22=abs(dFETy22)/max(abs(IFETy(2,2,:)));
% eFETy22=sum(adFETy22)/51;
% %误差形式2 max
                              dFETy11=ifety(1,1,:)-IFETy(1,1,:);
                              adFETy11=abs(dFETy11)/max(abs(IFETy(1,1,:)));
                              eFETy11=sum(adFETy11)/51;
                              dFETy12=ifety(1,2,:)-IFETy(1,2,:);
                              adFETy12=abs(dFETy12)/max(abs(IFETy(1,2,:)));
                              eFETy12=sum(adFETy12)/51;
                              dFETy21=ifety(2,1,:)-IFETy(2,1,:);
                              adFETy21=abs(dFETy21)/max(abs(IFETy(2,1,:)));
                              eFETy21=sum(adFETy21)/51;
                              dFETy22=ifety(2,2,:)-IFETy(2,2,:);
                              adFETy22=abs(dFETy22)/max(abs(IFETy(2,2,:)));
                              eFETy22=sum(adFETy22)/51;
%===
                              wy11=0.25;wy12=0.25;wy21=0.25;wy22=0.25;
                              FETyE=(wy11*eFETy11+wy12*eFETy12+wy21*eFETy21+wy22*eFETy22)*100;
% ===============ε1结束
  
%================计算总ε
                                      w1=0.5;             %2017年4月18日16:05:48， 只考虑S参数
                                      w2=0.5;
                                      e=w1*FETyE+w2*E;
%================结束

%                      %不同的判断标准， 记录robust提取的intrinsic参数值； 并记下误差，用于计算ε1+ε2  
%                      %SP为标准
%                          if E<RE
%                            RE=E;
%                            Var=var;
%                            FETyE;
%                            e11=es11;
%                            e12=es12;
%                            e21=es21;
%                            e22=es22;
%                            %2016年3月14日19:29:50  改进的waitbar  
% %                             waitbar(((i1-1)*n(2)*n(3)*n(4)*n(5)*n(6)*n(7)+(i2-1)*n(3)*n(4)*n(5)*n(6)*n(7)+(i3-1)*n(4)*n(5)*n(6)*n(7)+(i4-1)*n(5)*n(6)*n(7)+(i5-1)*n(6)*n(7)+(i6-1)*n(7)+i7)/(n(1)*n(2)*n(3)*n(4)*n(5)*n(6)*n(7)),h)
%                          end
                      %不同的判断标准， 记录robust提取的intrinsic参数值； 并记下误差，用于计算ε1+ε2  
                     %以总误差ε为标准 
                                      if e<Re
                                        Re=e;
                                        VarFETy=var;
                                        EfetS=E;
                                        yRE=FETyE;
                                          
                         
                                         
                                      end
                             end
                         end
                    end
               end
           end
       end
    end
end