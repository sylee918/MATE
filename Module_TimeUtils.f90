MODULE TIME_UTILS
   IMPLICIT NONE

CONTAINS

   pure elemental logical function is_leap_year(year)
      integer, intent(in) :: year
      is_leap_year = (mod(year, 4) == 0 .and. mod(year, 100) /= 0) .or. (mod(year, 400) == 0)
   end function is_leap_year

   pure elemental integer function days_in_year(year)
      integer, intent(in) :: year
      if (is_leap_year(year)) then
         days_in_year = 366
      else
         days_in_year = 365
      endif
   end function days_in_year

   ! Converts (year, doy) to continuous day count (Rata Die, days since 0001-01-01)
   pure elemental integer function ydoy_to_jday(year, doy)
      integer, intent(in) :: year, doy
      integer :: y
      y = year - 1
      ydoy_to_jday = y * 365 + y/4 - y/100 + y/400 + doy
   end function ydoy_to_jday

   ! Converts continuous day count back to (year, doy)
   pure subroutine jday_to_ydoy(jday, year, doy)
      integer, intent(in) :: jday
      integer, intent(out) :: year, doy
      integer :: y, start_jday, diy

      y = int(dble(jday - 1) / 365.2425d0) + 1
      start_jday = ydoy_to_jday(y, 1)

      do while (start_jday > jday)
         y = y - 1
         start_jday = ydoy_to_jday(y, 1)
      enddo

      diy = days_in_year(y)
      do while (jday >= start_jday + diy)
         y = y + 1
         start_jday = ydoy_to_jday(y, 1)
         diy = days_in_year(y)
      enddo

      year = y
      doy = jday - start_jday + 1
   end subroutine jday_to_ydoy

   ! Exact integer day difference between two YYYYDOYs (ydoy2 - ydoy1)
   pure elemental integer function ydoy_diff_days(ydoy2, ydoy1)
      integer, intent(in) :: ydoy2, ydoy1
      integer :: yr1, dy1, yr2, dy2
      yr1 = ydoy1 / 1000
      dy1 = mod(ydoy1, 1000)
      yr2 = ydoy2 / 1000
      dy2 = mod(ydoy2, 1000)
      ydoy_diff_days = ydoy_to_jday(yr2, dy2) - ydoy_to_jday(yr1, dy1)
   end function ydoy_diff_days

   ! Add integer days to integer YYYYDOY
   pure elemental integer function ydoy_add_days_int(ydoy, idays)
      integer, intent(in) :: ydoy, idays
      integer :: yr, dy, jday, new_jday, new_yr, new_dy
      yr = ydoy / 1000
      dy = mod(ydoy, 1000)
      jday = ydoy_to_jday(yr, dy)
      new_jday = jday + idays
      call jday_to_ydoy(new_jday, new_yr, new_dy)
      ydoy_add_days_int = new_yr * 1000 + new_dy
   end function ydoy_add_days_int

   ! Add real days to real YYYYDOY.fraction
   pure function ydoy_add_days(ydoy_real, ddays_real) result(new_ydoy_real)
      real*8, intent(in) :: ydoy_real, ddays_real
      real*8 :: new_ydoy_real
      integer :: yr, dy, jday, new_jday, new_yr, new_dy
      real*8 :: frac, total_days_frac
      integer :: int_shift

      yr = int(ydoy_real) / 1000
      dy = mod(int(ydoy_real), 1000)
      frac = ydoy_real - dble(int(ydoy_real))

      total_days_frac = frac + ddays_real
      int_shift = floor(total_days_frac)
      frac = total_days_frac - dble(int_shift)

      jday = ydoy_to_jday(yr, dy)
      new_jday = jday + int_shift
      call jday_to_ydoy(new_jday, new_yr, new_dy)
      new_ydoy_real = dble(new_yr * 1000 + new_dy) + frac
   end function ydoy_add_days

END MODULE TIME_UTILS
