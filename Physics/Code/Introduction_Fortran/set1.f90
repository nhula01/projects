implicit none

double precision :: F, C ! fahrenheit and celsius
double precision :: jan, feb, mar, apr, may, jun, jul, aug, sep, oct, nov, dec ! avg temp each month
double precision :: j_f_m_C,j_f_m_F, a_m_j_C,a_m_j_F,&
                    j_a_s_C,j_a_s_F, o_n_d_C,o_n_d_F  !avg temp each three months
                                !in celsisus and fahrenheit

!googling data we get
!the average data in fareinheit 
jan = 56.5 !68/45:High/Low consecutively
feb = 60.0
mar = 64.5
apr = 72.5
may = 80.5
jun = 90.0
jul = 94.0
aug = 93.0
sep = 86.5
oct = 76.5
nov = 63.5
dec = 56.0 

!calculating the average of each three months in fahrenheit
j_f_m_F = (jan+feb+mar)/3
a_m_j_F = (apr+may+jun)/3
j_a_s_F = (jul+aug+sep)/3
o_n_d_F = (oct+nov+dec)/3

! write the file
open(file="/Users/phihung/Documents/PHY/first/homework/hw_01/temp.dat", unit=10)
    ! call subroutine to do the conversions
    call sub_fun(j_f_m_C,j_f_m_F)
    call sub_fun(a_m_j_C,a_m_j_F) 
    call sub_fun(j_a_s_C,j_a_s_F) 
    call sub_fun(o_n_d_C,o_n_d_F)  
    ! write to the file the data
    write(10,*) j_f_m_C,j_f_m_F
    write(10,*) a_m_j_C,a_m_j_F
    write(10,*) j_a_s_C,j_a_s_F
    write(10,*) o_n_d_C,o_n_d_F

close(10)
end

! subroutine for conversions
subroutine sub_fun(C,F)
    implicit none
    double precision C
    double precision F 
    C = (F-32)/1.8
end subroutine sub_fun

