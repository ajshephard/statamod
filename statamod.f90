
!    This file is part of STATAMOD. Copyright (c) 2006-2014 Andrew Shephard
!
!    STATAMOD is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    STATAMOD is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with STATAMOD.  If not, see <http://www.gnu.org/licenses/>.

module stataMod

    implicit none

    private

    !public access
    public :: openStata, saveStata, closeOpenStata, closeSaveStata, readStata
    public :: writeStata, descStata, nobsStata, nvarStata, existStata
    public :: maxST_byte, maxST_int, maxST_long, maxST_float, maxST_double
    public :: minST_byte, minST_int, minST_long, minST_float, minST_double
    public :: sampleStata

    !data types
    integer, parameter :: sp = selected_real_kind(6,30)
    integer, parameter :: dp = selected_real_kind(15,100)
    integer, parameter :: i4 = selected_int_kind(5)
    integer, parameter :: i2 = selected_int_kind(3)
    integer, parameter :: i1 = selected_int_kind(1)

    integer,  parameter :: maxST_byte   = 100
    integer,  parameter :: maxST_int    = 32740
    integer,  parameter :: maxST_long   = 2147483620
    real(sp), parameter :: maxST_float  = Z'7effffff'         !approx +1.701e+38
    real(dp), parameter :: maxST_double = Z'7fdfffffffffffff' !approx +8.988e+307

    integer,  parameter :: minST_byte   = -127
    integer,  parameter :: minST_int    = -32767
    integer,  parameter :: minST_long   = -2147483647
    real(sp), parameter :: minST_float  = Z'feffffff'         !approx -1.701e+38
    real(dp), parameter :: minST_double = Z'ffefffffffffffff' !approx -1.798e+308

    type :: position
        integer :: start
        integer, dimension(:), allocatable :: offset
        integer, dimension(:), allocatable :: varpos
        integer :: nvaroffset
    end type position

    type stcell !use linked list when saving
        real(dp),    dimension(:), allocatable :: valReal8
        real(sp),    dimension(:), allocatable :: valReal4
        integer(i4), dimension(:), allocatable :: valInt4
        integer(i2), dimension(:), allocatable :: valInt2
        integer(i1), dimension(:), allocatable :: valInt1
        integer(i2)                            :: valType
        character(33)                          :: name
        character(81)                          :: label
        type(stcell), pointer                  :: next
    end type stcell

    type :: stata

        ! file header
        integer(i1)   :: ds_format
        integer(i1)   :: byteorder
        integer(i1)   :: filetype
        integer(i2)   :: nvar
        integer(i4)   :: nobs
        character(81) :: data_label
        character(18) :: time_stamp

        ! descriptors
        integer(i2),   dimension(:), allocatable :: typlist
        character(33), dimension(:), allocatable :: varlist
        integer(i2),   dimension(:), allocatable :: srtlist
        character(49), dimension(:), allocatable :: fmtlist
        character(33), dimension(:), allocatable :: lbllist

        !variable labels
        character(81), dimension(:),   allocatable :: varlabl
        character(1),  dimension(:,:), allocatable :: thesedata !enitre data if cached

        integer :: unit  !unit number
        logical :: cache !cache data?
        logical :: init = .false. !initialised?

        logical :: doSample = .false.
        integer, dimension(:), allocatable :: sampleWeight

        type(position) :: stataPosition  !file position

        ! the remaining information is for when Stata files are written

        integer(i4) :: saveDim   !array dim size
        integer(i4) :: saveUnit  !unit number
        integer(i4) :: saveNObs  !save observations
        integer(i2) :: saveNVar  !save variables

        logical :: saveInit = .false. !initialised?
        logical :: saveOnce  !has a variable been written?
        logical :: saveCache !save cache?
        logical :: saveCompress !compress dataset

        character(81) :: saveLabel !dataset label
        character(18) :: saveTime  !time stamp

        type(stcell), pointer :: saveCurr, saveHead

    end type stata

    type(stata), save :: statafile

    interface writestata
        module procedure writeStataInt1
        module procedure writeStataInt2
        module procedure writeStataInt4
        module procedure writeStataReal4
        module procedure writeStataReal8
    end interface writestata

    interface readStata
        module procedure readStataInt1
        module procedure readStataInt2
        module procedure readStataInt4
        module procedure readStataReal4
        module procedure readStataReal8
        module procedure readStataLogical
    end interface readStata

contains

    subroutine sampleStata(sampleWeight)
        implicit none
        integer(i4), intent(in), optional :: sampleWeight(:)

        if (.not. stataFile%init) then
            call statamodError('STATA file is not open')
        end if

        if (present(sampleWeight)) then
            if (size(sampleWeight) .ne. stataFile%nobs) then
                call statamodError('sampleWeight dimension does not equal number of observations')
            end if
            if ( (minval(sampleWeight)<1) .or. (maxval(sampleWeight)>stataFile%nobs) ) then
                call statamodError('sampleWeight is out of range')
            end if
            stataFile%doSample = .true.
            if (allocated(stataFile%sampleWeight)) then
                deallocate(stataFile%sampleWeight)
            end if

            allocate(stataFile%sampleWeight(stataFile%nobs))

            stataFile%sampleWeight = sampleWeight
        else
            stataFile%doSample = .false.
        end if

    end subroutine sampleStata

    subroutine saveStata(fileName,obs,label)

        implicit none

        character(*),  intent(in)           :: fileName
        integer(i4),   intent(in)           :: obs
        character(*),  intent(in), optional :: label

        integer(i1) :: tempInt1
        integer(i4) :: ios

        if (stataFile%saveInit) then
            call statamodError('A file has already been specified for saving')
        end if

        if (obs < 1) then
            call statamodError('Expecting at least one observation to be declared')
        end if

        stataFile%saveNObs = obs
        stataFile%saveNVar = 0
        stataFile%saveOnce = .false.
        stataFile%saveTime = dateStata()

        if (present(label)) then
            tempInt1 = len(label) + 1
            stataFile%saveLabel = label
            stataFile%saveLabel(tempInt1:81) = char(0)
        else
            stataFile%saveLabel = char(0)
        end if

        call getUnit(stataFile%saveUnit)
        open(UNIT=stataFile%saveUnit,FILE=fileName,ACTION='write',STATUS='REPLACE',ACCESS='stream',IOSTAT=ios)

        if (ios /= 0) then
            call statamodError('error saving file '//fileName)
        end if

        !save as stata 7 s/e
        write(stataFile%saveUnit, IOSTAT=ios) 111_i1 ! Stata 7 S/E
        write(stataFile%saveUnit, IOSTAT=ios) 2_i1   ! byte order
        write(stataFile%saveUnit, IOSTAT=ios) 1_i1   ! filetype
        write(stataFile%saveUnit, IOSTAT=ios) 0_i1   ! junk

        if (ios /= 0) then
            call statamodError('error writing file header')
        end if

        stataFile%saveInit = .true.
        stataFile%saveDim  = 0

    end subroutine saveStata

    subroutine writeStataReal4(thisVar, thisName, thisLabel)

        implicit none

        real(sp), dimension(:), intent(in) :: thisVar

        include 'writestata1.inc'

        ! var type
        stataFile%saveCurr%valType = 254
        stataFile%saveDim = stataFile%saveDim + 4

        ! copy contents
        allocate(stataFile%saveCurr%valReal4(stataFile%saveNObs), STAT=ios)

        if (ios /= 0) then
            call statamodError('error allocating memory')
        end if

        stataFile%saveCurr%valReal4 = thisVar

        include 'writestata2.inc'

    end subroutine writeStataReal4

    subroutine writeStataReal8(thisVar, thisName, thisLabel)

        implicit none

        real(dp), dimension(:), intent(in) :: thisVar

        include 'writestata1.inc'

        ! var type
        stataFile%saveCurr%valType = 255
        stataFile%saveDim = stataFile%saveDim + 8

        ! copy contents
        allocate(stataFile%saveCurr%valReal8(stataFile%saveNObs), STAT=ios)

        if (ios /= 0) then
            call statamodError('error allocating memory')
        end if

        stataFile%saveCurr%valReal8 = thisVar

        include 'writestata2.inc'

    end subroutine writeStataReal8

    subroutine writeStataInt1(thisVar, thisName, thisLabel)

        implicit none

        integer(i1), dimension(:), intent(in) :: thisVar

        include 'writestata1.inc'

        ! var type
        stataFile%saveCurr%valType = 251
        stataFile%saveDim = stataFile%saveDim + 1

        ! copy contents
        allocate(stataFile%saveCurr%valInt1(stataFile%saveNObs), STAT=ios)

        if (ios /= 0) then
            call statamodError('error allocating memory')
        end if

        stataFile%saveCurr%valInt1 = thisVar

        include 'writestata2.inc'

    end subroutine writeStataInt1

    subroutine writeStataInt2(thisVar, thisName, thisLabel)

        implicit none

        integer(i2), dimension(:), intent(in) :: thisVar

        include 'writestata1.inc'

        ! var type
        stataFile%saveCurr%valType = 252
        stataFile%saveDim = stataFile%saveDim + 2

        ! copy contents
        allocate(stataFile%saveCurr%valInt2(stataFile%saveNObs), STAT=ios)

        if (ios /= 0) then
            call statamodError('error allocating memory')
        end if

        stataFile%saveCurr%valInt2 = thisVar

        include 'writestata2.inc'

    end subroutine writeStataInt2

    subroutine writeStataInt4(thisVar, thisName, thisLabel)

        implicit none

        integer(i4), dimension(:), intent(in) :: thisVar

        include 'writestata1.inc'

        ! var type
        stataFile%saveCurr%valType = 253
        stataFile%saveDim = stataFile%saveDim + 4

        ! copy contents
        allocate(stataFile%saveCurr%valInt4(stataFile%saveNObs), STAT=ios)

        if (ios /= 0) then
            call statamodError('error allocating memory')
        end if

        stataFile%saveCurr%valInt4 = thisVar

        include 'writestata2.inc'

    end subroutine writeStataInt4

    subroutine closeSaveStata(cache)

        implicit none

        logical, optional :: cache

        integer(i4)  :: ios, tempInt4, i, j
        character(1) :: unsignedByte !unsigned integers are not supported. use this as work around

        integer(i1),   dimension(2*(stataFile%saveNVar+1)) :: srtlist
        character(12), dimension(stataFile%saveNVar)       :: fmtlist
        character(33), dimension(stataFile%saveNVar)       :: lbllist

        integer(i1), dimension(5)   :: expfield

        character(stataFile%saveDim), allocatable, dimension(:) :: dataset

        if (.not. stataFile%saveInit) then
            call statamodError('no file have been specified for saving')
        end if

        write(stataFile%saveUnit, IOSTAT=ios) stataFile%saveNVar  !Variables
        write(stataFile%saveUnit, IOSTAT=ios) stataFile%saveNObs  !Observations
        write(stataFile%saveUnit, IOSTAT=ios) stataFile%saveLabel !Label
        write(stataFile%saveUnit, IOSTAT=ios) stataFile%saveTime  !time stamp

        if (ios /= 0) then
            call statamodError('error writing file specific header information')
        end if

        ! write descriptors

        ! typlist
        stataFile%saveCurr => stataFile%saveHead
        tempInt4 = 0

        do i = 1, stataFile%saveNVar
            unsignedByte = char(stataFile%saveCurr%valType) !work around
            write(stataFile%saveUnit, IOSTAT=ios) unsignedByte
            stataFile%saveCurr => stataFile%saveCurr%next
        end do !i

        ! varlist
        stataFile%saveCurr => stataFile%saveHead

        do i = 1, stataFile%saveNVar
            write(stataFile%saveUnit, IOSTAT=ios) stataFile%saveCurr%name
            stataFile%saveCurr => stataFile%saveCurr%next
        end do !i

        ! srtlist (unsorted)
        srtlist = 0
        write(stataFile%saveUnit, IOSTAT=ios) srtlist

        ! fmtlist (use default stata format)
        fmtlist = '%8.0g'//char(0)
        write(stataFile%saveUnit, IOSTAT=ios) fmtlist

        ! lbllist
        lbllist = char(0)
        write(stataFile%saveUnit, IOSTAT=ios) lbllist

        if (ios /= 0) then
            call statamodError('error writing descriptors')
        end if

        ! write variable labels

        stataFile%saveCurr => stataFile%saveHead

        do i = 1, stataFile%saveNVar
            write(stataFile%saveUnit, IOSTAT=ios) stataFile%saveCurr%label
            stataFile%saveCurr => stataFile%saveCurr%next
        end do !i

        if (ios /= 0) then
            call statamodError('error writing variable labels')
        end if

        ! write expansion fields

        expfield = 0
        write(stataFile%saveUnit, IOSTAT=ios) expfield

        if (ios /= 0) then
            call statamodError('error writing expansion fields')
        end if

        ! write actual variables

        if (present(cache)) then
            stataFile%saveCache = cache
        else
            stataFile%saveCache = .true.
        end if

        if (stataFile%saveCache) then

            ! because the data is of potentially different types, we must construct
            ! an appropriately sized character array and use the transfer function
            ! to copy exact bit pattern

            allocate(dataset(stataFile%saveNObs), STAT=ios)

            if (ios /= 0) then
                call statamodError('error allocating memory for variable writing. Use closeSaveStata(.false.).')
            end if

            stataFile%saveCurr => stataFile%saveHead
            tempInt4 = 1

            do i = 1, stataFile%saveNVar

                select case(stataFile%saveCurr%valType)
                    case(251_i2)   !ST_byte
                        dataset(:)(tempInt4:tempInt4)   = transfer(stataFile%saveCurr%valInt1,(/1_'0'/))
                        tempInt4 = tempInt4 + 1
                    case(252_i2)   !ST_int
                        dataset(:)(tempInt4:tempInt4+1) = transfer(stataFile%saveCurr%valInt2,(/1_'00'/))
                        tempInt4 = tempInt4 + 2
                    case(253_i2)   !ST_long
                        dataset(:)(tempInt4:tempInt4+3) = transfer(stataFile%saveCurr%valInt4,(/1_'0000'/))
                        tempInt4 = tempInt4 + 4
                    case(254_i2)   !ST_float
                        dataset(:)(tempInt4:tempInt4+3) = transfer(stataFile%saveCurr%valReal4,(/1_'0000'/))
                        tempInt4 = tempInt4 + 4
                    case(255_i2)   !ST_double
                        dataset(:)(tempInt4:tempInt4+7) = transfer(stataFile%saveCurr%valReal8,(/1_'00000000'/))
                        tempInt4 = tempInt4 + 8
                    case default
                        call statamodError('unknown data type')
                end select

                stataFile%saveCurr => stataFile%saveCurr%next

            end do !i

            write(stataFile%saveUnit, IOSTAT=ios) dataset

        else

            do i = 1, stataFile%saveNObs

                stataFile%saveCurr => stataFile%saveHead

                do j = 1, stataFile%saveNVar

                    select case(stataFile%saveCurr%valType)
                        case(251_i2)   !ST_byte
                            write(stataFile%saveUnit,IOSTAT=ios) stataFile%saveCurr%valInt1(i:i)
                        case(252_i2)   !ST_int
                            write(stataFile%saveUnit,IOSTAT=ios) stataFile%saveCurr%valInt2(i:i)
                        case(253_i2)   !ST_long
                            write(stataFile%saveUnit,IOSTAT=ios) stataFile%saveCurr%valInt4(i:i)
                        case(254_i2)   !ST_float
                            write(stataFile%saveUnit,IOSTAT=ios) stataFile%saveCurr%valReal4(i:i)
                        case(255_i2)   !ST_double
                            write(stataFile%saveUnit,IOSTAT=ios) stataFile%saveCurr%valReal8(i:i)
                        case default
                            call statamodError('unknown data type')
                    end select

                stataFile%saveCurr => stataFile%saveCurr%next

                end do !j

            end do !i

        end if

        if (ios /= 0) then
            call statamodError('error writing actual data contents')
        end if

        ! deallocate memory and close file

        stataFile%saveCurr => stataFile%saveHead

        do while(associated(stataFile%saveCurr))
            stataFile%saveHead => stataFile%saveCurr%next
            deallocate(stataFile%saveCurr, STAT=ios)
            if (ios /= 0) call statamodWarn('error deallocating memory')
            stataFile%saveCurr => stataFile%saveHead
        end do

        stataFile%saveInit = .false.

        !if (stataFile%saveCache) then
        !    deallocate(dataset)
        !end if

        !if (allocated(dataset)) deallocate(dataset)

        close(stataFile%saveUnit)

        call statamodMsg('STATA file successfully saved')

    end subroutine closeSaveStata


    subroutine openStata(fileName,cache)

        implicit none

        character(*), intent(in)           :: fileName
        logical,      intent(in), optional :: cache

        logical       :: ex
        integer(i4)   :: ios
        integer(i4)   :: i, j

        integer(i1)   :: tempInt1 !temporary integers
        integer(i4)   :: tempInt4

        character(9)  :: str9 !temporary strings
        character(32) :: str32
        character(33) :: str33
        character(81) :: str81

        logical      :: cacheData = .false.
        character(1) :: unsignedByte !unsigned integers are not supported. use this as work around

        integer :: filePos ! file position

        if (stataFile%init) then
            call statamodError('STATA file already open')
        end if

        stataFile%init = .false. !change to true if successful

        ! open file
        inquire(FILE=fileName,EXIST=ex)
        if (.not. ex) then
            call statamodError('file '//fileName//' does not exist')
        end if

        call getUnit(stataFile%unit)
        open(UNIT=stataFile%unit,FILE=fileName,action='read',status='old',access='stream',iostat=ios)

        ! error check
        if (ios/=0) then
            call statamodError('error opening file '//filename)
        end if

        ! read header
        read(stataFile%unit, IOSTAT=ios) stataFile%ds_format
        read(stataFile%unit, IOSTAT=ios) stataFile%byteorder
        read(stataFile%unit, IOSTAT=ios) stataFile%filetype
        read(stataFile%unit, IOSTAT=ios) tempInt1 !junk
        read(stataFile%unit, IOSTAT=ios) stataFile%nvar
        read(stataFile%unit, IOSTAT=ios) stataFile%nobs

        if (versionStata(stataFile%ds_format)=='5') then
            read(stataFile%unit, IOSTAT=ios) str32
            stataFile%data_label = str32
        else
            read(stataFile%unit, IOSTAT=ios) str81
            stataFile%data_label = str81
        end if

        read(stataFile%unit, IOSTAT=ios) stataFile%time_stamp

        if (ios==1) then
            call statamodError('error reading from file '//fileName)
        end if

        !check header
        if (versionStata(stataFile%ds_format)=='?') ios=1
        if (stataFile%byteorder/=1 .and. stataFile%byteorder/=2) ios=1
        if (stataFile%filetype/=1) ios=1
        if (stataFile%nvar<0) ios=1
        if (stataFile%nobs<0) ios=1

        if (ios==1) then
            call statamodError('not a valid STATA file')
        end if

        if (stataFile%nvar==0) then
            call statamodError('file does not contain any variables')
        end if

        if (stataFile%nobs==0) then
            call statamodError('file does not contain any observations')
        end if

        !read descriptors
        allocate(stataFile%typlist(stataFile%nvar), STAT=ios)
        allocate(stataFile%varlist(stataFile%nvar), STAT=ios)
        allocate(stataFile%srtlist(stataFile%nvar + 1), STAT=ios)
        allocate(stataFile%fmtlist(stataFile%nvar), STAT=ios)
        allocate(stataFile%lbllist(stataFile%nvar), STAT=ios)

        if (ios /= 0) then
            call statamodError('error allocating memory')
        end if

        do i = 1, stataFile%nvar
            read(stataFile%unit) unsignedByte
            stataFile%typlist(I) = ichar(unsignedByte)
        end do

        do i = 1, stataFile%nvar

            if (versionStata(stataFile%ds_format)=='5' .or. versionStata(stataFile%ds_format)=='6') then
                read(stataFile%unit) str9
                stataFile%varlist(i) = str9
            else
                read(stataFile%unit) str33
                stataFile%varlist(i) = str33
            end if

            do j = 1, 33
                if (stataFile%varlist(i)(j:j) == char(0)) then
                        stataFile%varlist(i)(j:33) = ' '
                        exit
                end if
            end do

        end do

        do i = 1, (stataFile%nvar + 1)
            read(stataFile%unit) stataFile%srtlist(i)
        end do

        if (stataFile%ds_format>=114) then
            do i = 1, stataFile%nvar
                read(stataFile%unit) stataFile%fmtlist(i)
            end do
        else
            do i = 1, stataFile%nvar
                read(stataFile%unit) stataFile%fmtlist(i)(1:12)
            end do
        end if

        do i = 1, stataFile%nvar

            read(stataFile%unit) stataFile%lbllist(i)
        end do

        !read variable labels

        allocate(stataFile%varlabl(stataFile%nvar), STAT=ios)

        if (ios /= 0) then
            call statamodError('error allocating memory')
        end if

        do i = 1, stataFile%nvar
            read(stataFile%unit) stataFile%varlabl(i);
        end do

        !skip expansion fields

        tempInt1 = 1
        tempInt4 = 1

        inquire(unit=stataFile%unit, POS=filePos)

        do while ( tempInt1 /= 0 .and. tempInt4 /=0 )
            read(stataFile%unit,POS=filePos) tempInt1
            read(stataFile%unit) tempInt4
            inquire(unit=stataFile%unit, POS=filePos)
            filePos = filePos + tempInt4
        end do
        !position of start of data

        inquire(unit=stataFile%unit, POS=stataFile%stataPosition%start)

        allocate(stataFile%stataPosition%offset(stataFile%nvar))
        allocate(stataFile%stataPosition%varpos(stataFile%nvar))

        stataFile%stataPosition%offset(1) = 0
        stataFile%stataPosition%varpos(1) = stataFile%stataPosition%start

        if (stataFile%nvar > 1) then
            do i = 2, stataFile%nvar
                stataFile%stataPosition%offset(i) = stataFile%stataPosition%offset(i-1) + sizeStata(stataFile%typlist(i-1))
                stataFile%stataPosition%varpos(i) = stataFile%stataPosition%varpos(i-1) + sizeStata(stataFile%typlist(i-1))
            end do
        end if

        stataFile%stataPosition%nvarOffset = stataFile%stataPosition%offset(stataFile%nvar) + &
            & sizeStata(stataFile%typlist(stataFile%nvar))

        !load entire dataset into memory

        if (present(cache)) then
            if (cache) cacheData = .true.
        else
            cacheData = .true.
        end if

        if (cacheData) then

            allocate(stataFile%theseData(stataFile%stataPosition%nvarOffset,stataFile%nobs), STAT = ios)

            if (ios /= 0) then
                call statamodError('error allocating memory')
            end if

            read(stataFile%unit, IOSTAT = ios) stataFile%theseData

            if (ios /=0) then
                call statamodError('error reading data into memory')
            end if

            stataFile%cache = .true.
            close(stataFile%unit)

        end if

        !if we reach here, data has opened okay
        stataFile%init = .true.

    end subroutine openStata

    subroutine readStataReal4(readStataVar,varname)

        implicit none

        character(*), intent(in)                :: varname
        real(sp),     intent(out), dimension(:) :: readStataVar

        real(sp),  dimension(size(readStataVar)) :: readStataVar2

        include 'readstata.inc'

    end subroutine readStataReal4

    subroutine readStataReal8(readStataVar,varname)

        implicit none

        character(*), intent(in)                :: varname
        real(dp),     intent(out), dimension(:) :: readStataVar

        real(dp),  dimension(size(readStataVar)) :: readStataVar2

        include 'readstata.inc'

    end subroutine readStataReal8

    subroutine readStataInt1(readStataVar,varname)

        implicit none

        character(*), intent(in)                :: varname
        integer(i1),  intent(out), dimension(:) :: readStataVar

        integer(i1),  dimension(size(readStataVar)) :: readStataVar2

        include 'readstata.inc'

    end subroutine readStataInt1

    subroutine readStataInt2(readStataVar,varname)

        implicit none

        character(*), intent(in)                :: varname
        integer(i2),  intent(out), dimension(:) :: readStataVar

        integer(i2),  dimension(size(readStataVar)) :: readStataVar2

        include 'readstata.inc'

    end subroutine readStataInt2

    subroutine readStataInt4(readStataVar,varname)

        implicit none

        character(*), intent(in)                :: varname
        integer(i4),  intent(out), dimension(:) :: readStataVar

        integer(i4),  dimension(size(readStataVar)) :: readStataVar2

        include 'readstata.inc'

    end subroutine readStataInt4

    subroutine readStataLogical(readStataVarLogical,varname)

        implicit none

        character(*), intent(in)                :: varname
        logical,      intent(out), dimension(:) :: readStataVarLogical
        real(dp),     dimension(size(readStataVarLogical)) :: readStataVar
        real(dp),     dimension(size(readStataVarLogical)) :: readStataVar2

        include 'readstata.inc'

        where (readStataVar>0.0_dp)
            readStataVarLogical = .true.
        else where
            readStataVarLogical = .false.
        end where

    end subroutine readStataLogical

    subroutine closeOpenStata()

        implicit none

        integer(i4) :: err

        err = 0

        if (.not. stataFile%init) then
            call statamodWarn('STATA file is not open')
            return
        end if

        if (allocated(stataFile%theseData)) deallocate(stataFile%theseData, STAT=err)

        if (allocated(stataFile%typlist))   deallocate(stataFile%typlist,   STAT=err)
        if (allocated(stataFile%varlist))   deallocate(stataFile%varlist,   STAT=err)
        if (allocated(stataFile%srtlist))   deallocate(stataFile%srtlist,   STAT=err)
        if (allocated(stataFile%fmtlist))   deallocate(stataFile%fmtlist,   STAT=err)
        if (allocated(stataFile%lbllist))   deallocate(stataFile%lbllist,   STAT=err)
        if (allocated(stataFile%varlabl))   deallocate(stataFile%varlabl,   STAT=err)

        if (allocated(stataFile%stataPosition%offset)) deallocate(stataFile%stataPosition%offset, STAT=err)

        if (err /= 0) call statamodWarn('error deallocating memory')

        if (.not. stataFile%cache) close(stataFile%unit) !already closed if cache is true

        stataFile%init = .false.

    end subroutine closeOpenStata

    function existStata(varname) !check whether var exists in dataset

        implicit none

        character(*), intent(in) :: varname
        logical                  :: existStata

        integer(i4) :: thisVarno

        if (.not. stataFile%init) then
            call statamodWarn('STATA file is not open')
            existStata = .false.
            return
        end if

        thisVarno = varno(varname)

        if (thisVarno==0) then
            existStata = .false.
        else
            existStata = .true.
        end if

    end function existStata

    integer(i4) function nobsStata() !return number of observations

        implicit none

        if (.not. stataFile%init) then
            call statamodWarn('STATA file is not open')
            nobsStata = 0
            return
        end if

        nobsStata = stataFile%nobs

    end function nobsStata


    integer(i2) function nvarStata() !return number of variables

        implicit none

        if (.not. stataFile%init) then
            call statamodWarn('STATA file is not open')
            nvarStata = 0
            return
        end if

        nvarStata = stataFile%nvar

    end function nvarStata


    subroutine descStata() !print descriptive statistics to console

        implicit none

        character(24) :: intToStr
        integer(i4)   :: i

        if (.not. stataFile%init) then
            call statamodWarn('STATA file is not open')
            return
        end if

        write (*,*)  'ST_version: ' // versionStata(stataFile%ds_format)
        write (*,*)  'ST_date:    ' // stataFile%time_stamp

        write (intToStr,*) stataFile%nvar
        write (*,*)  'ST_nvars:   ' // adjustl(intToStr)

        write (intToStr,*) stataFile%nobs
        write (*,*)  'ST_nobs:    ' // adjustl(intToStr)

        write (*,*)

        do i = 1, stataFile%nvar
            write (*,*)  typeStata(i) // stataFile%varlist(i)
        end do

    end subroutine descStata


    character(12) function typeStata(varno) !return string with typlist description

        implicit none

        integer(i4), intent(in) :: varno

        select case(stataFile%typlist(varno))
            case(1_i2:244_i2) !str
                typeStata = 'ST_str'
            case(251_i2)   !ST_byte
                typeStata = 'ST_byte'
            case(252_i2)   !ST_int
                typeStata = 'ST_int'
            case(253_i2)   !ST_long
                typeStata = 'ST_long'
            case(254_i2)   !ST_float
                typeStata = 'ST_float'
            case(255_i2)   !ST_double
                typeStata = 'ST_double'
            case default
                typeStata = 'unknown'
        end select

    end function typeStata


    integer(i2) function sizeStata(typlist) !return size of typlist

        implicit none

        integer(i2), intent(in) :: typlist

        select case(typlist)
            case(1_i2:244_i2) !str
                sizeStata = typlist
            case(251_i2)   !ST_byte
                sizeStata = 1_i2
            case(252_i2)   !ST_int
                sizeStata = 2_i2
            case(253_i2)   !ST_long
                sizeStata = 4_i2
            case(254_i2)   !ST_float
                sizeStata = 4_i2
            case(255_i2)   !ST_double
                sizeStata = 8_i2
            case default
                sizeStata = 0_i2
        end select

    end function sizeStata

    character(3) function versionStata(version) !return string with STATA version number

        implicit none

        integer(i1), intent(in) :: version

        select case(version)
            case(114_i1)
                versionStata = '10'
            case(113_i1)
                versionStata = '8'
            case(111_i1)
                versionStata = '7SE'
            case(110_i1)
                versionStata = '7'
            case(108_i1)
                versionStata = '6'
            case(105_i1)
                versionStata = '5'
            case default
                versionStata = '?'
        end select

    end function versionStata

    integer(i2) function varno(varname) !return variable number corresponding to varname

        implicit none

        character(*), intent(in) :: varname
        integer(i4)              :: i

        varno = 0

        do i = 1, stataFile%nvar
            if (varname == stataFile%varlist(i)) then
                varno = i
                exit
            end if
        end do

    end function varno

    character(18) function dateStata()

        implicit none

        character(18) :: date_time(2)
        character(3)  :: this_month

        call date_and_time(date_time(1), date_time(2))

        select case(date_time(1)(5:6))
            case('01')
                this_month = 'Jan'
            case('02')
                this_month = 'Feb'
            case('03')
                this_month = 'Mar'
            case('04')
                this_month = 'Apr'
            case('05')
                this_month = 'May'
            case('06')
                this_month = 'Jun'
            case('07')
                this_month = 'Jul'
            case('08')
                this_month = 'Aug'
            case('09')
                this_month = 'Sep'
            case('10')
                this_month = 'Oct'
            case('11')
                this_month = 'Nov'
            case('12')
                this_month = 'Dec'
            case default
                this_month = '???'
        end select

        dateStata = date_time(1)(7:8)//' '//this_month//' '//date_time(1)(1:4)//' '//date_time(2)(1:2)//':'//date_time(2)(3:4)
        dateStata(18:18) = char(0)

    end function dateStata

    subroutine statamodError(errmsg,funit)

        implicit none

        character(*),      intent(in) :: errmsg
        integer, optional, intent(in) :: funit

        if (present(funit)) then
            write (funit,*) 'STATAMOD ERROR: ',errmsg
            stop 'program terminated by statamodError'
        else
            write (*,*) 'STATAMOD ERROR: ',errmsg
            stop 'program terminated by statamodError'
        end if

        return

    end subroutine statamodError

    subroutine statamodWarn(warnmsg,funit)

        character(*),      intent(in) :: warnmsg
        integer, optional, intent(in) :: funit

        if (present(funit)) then
            write (funit,*) 'STATAMOD WARNING: ', warnmsg
        else
            write (*,*) 'STATAMOD WARNING: ', warnmsg
        end if

        return

    end subroutine statamodWarn

    subroutine statamodMsg(msg,funit)

        character(*),      intent(in) :: msg
        integer, optional, intent(in) :: funit

        if (present(funit)) then
            write (funit,*) 'STATAMOD: ', msg
        else
            write (*,*) 'STATAMOD: ', msg
        end if

        return

    end subroutine statamodMsg

    subroutine getUnit(funit)
        use, intrinsic :: iso_fortran_env
        implicit none
        integer, intent(out) :: funit
        integer        :: i
        logical        :: opend
        integer, parameter   :: stdout = output_unit
        integer, parameter   :: maxunit = 99
        do i = stdout+1, maxunit
            inquire(unit=i, opened=opend)
            if (.not. opend) then
                funit = i
                exit
            endif
        end do
    end subroutine getUnit

end module stataMod
