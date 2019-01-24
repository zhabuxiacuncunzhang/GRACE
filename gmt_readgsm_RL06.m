function    [cs,cs_sigma,int_year,int_month,time] = gmt_readgsm_RL06(dir_in,file_name,lmax)

% Read the gravity filed files
%
% INPUT:
%   dir_in      full path
%   file_name   filename 
%   lmax        maximum degree in files
% 
% OUTPUT:
%   cs          spherical harmonic coefficients in CS-format
%   cs_sigma    spherical harmonic coefficients in CS-format (formal error)
%   int_year    year
%   int_month   mont
%   time        fraction in the year
%
cs=zeros(lmax+1,lmax+1);
cs_sigma=zeros(lmax+1,lmax+1);

 % read the header
    head_index=0;
    fid = fopen([dir_in file_name],'r');
    tline = fgetl(fid);
    head_index=125;
    degree_max=lmax;
    fclose(fid);
  [~, l, m, Clm, Slm, Clm_sigma, Slm_sigma] = textread(strcat(dir_in,file_name),'%s %u %u %f %f %f %f %*[^\n]','headerlines',head_index);
    for i = 1:length(l)
        sc_tmp(l(i)+1,degree_max+1-m(i)) = Slm(i);
        sc_tmp(l(i)+1,degree_max+1+m(i)) = Clm(i);
        sc_sigma_tmp(l(i)+1,degree_max+1-m(i)) = Slm_sigma(i);
        sc_sigma_tmp(l(i)+1,degree_max+1+m(i)) = Clm_sigma(i);
    end
    cs_tmp=gmt_sc2cs(sc_tmp);
    cs_sigma_tmp=gmt_sc2cs(sc_sigma_tmp);
    % Get SH coefficients
    cs = cs_tmp(1:lmax+1,1:lmax+1);
    cs_sigma= cs_sigma_tmp(1:lmax+1,1:lmax+1);
    % Get time tag GSM-2_2003001-2003031_0026_EIGEN_G---_0005
    year1 = str2num(file_name(7:10));
    year2 = str2num(file_name(15:18));
    day1  = str2num(file_name(11:13));
    day2  = str2num(file_name(19:21));
    if year1 == year2
        meanday = (day1+day2)/2;
    else
        if (day1+(366-day1+day2)/2)>365 % in latter year
            year1   = year1 + 1;
            meanday = day2-(366-day1+day2)/2;
        else
            meanday = day1+(366-day1+day2)/2;  % in former year
        end
    end
    time        = year1 + meanday/365.;
    meanday     = round(meanday);
    int_year    = year1;
    int_month   = gmt_get_mon_day(year1,meanday+1);
end
