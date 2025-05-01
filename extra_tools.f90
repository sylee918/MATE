



      Subroutine Generate_tag(lon,lat, tag)
         ! Generate "tag" in format "i3.3" considering negative latitudes.
         ! Example: tag = "_lon270_lat000"
         USE SETTING
         real*8 lon, lat
         character*30 tag

         if (lat .ge. 0) then
            write(tag,'(A, i3.3, A, i3.3, A)') "_lon", nint(lon*180/pi), "_lat", nint(lat*180/pi)
         else
            write(tag,'(A, i3.3, A, i3.2, A)') "_lon", nint(lon*180/pi), "_lat", nint(lat*180/pi)
         endif

         return
      End







