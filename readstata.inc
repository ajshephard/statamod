
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

integer(i1) :: tempST_byte
integer(i2) :: tempST_int
integer(i4) :: tempST_long
real(sp)    :: tempST_float
real(dp)    :: tempST_double

integer(i4) :: i, thisVarno, varStart, sizeReadStata
integer     :: filePos


if (.not. stataFile%init) then
    call statamodError('STATA file is not open')
end if

!remove whitespace and find index
thisVarno = varno(trim(adjustl(varname)))

if (thisVarno==0) then
    call statamodWarn('variable '//trim(adjustl(varname))//' does not exist')
    readStataVar = 0
    return
end if

sizeReadStata = size(readStataVar)

if (sizeReadStata.lt.stataFile%nobs) then
    call statamodWarn('number of observations exceeds array storage')
end if

!for access from memory

if (stataFile%cache) then

    varStart = stataFile%stataPosition%offset(thisVarno) + 1 ! column location of varname

    !if (sizeReadStata==stataFile%nobs) then
        !use transfer intrinsic to cast data
        select case(stataFile%typlist(thisVarno))
            case(1_i2:244_i2) !str
                call statamodWarn('strings are not supported; setting array to zero')
                readStataVar(1:stataFile%nobs) = 0
            case(251_i2)   !ST_byte
                readStataVar(1:stataFile%nobs) = transfer(stataFile%theseData(varStart:varStart,:),(/0_i1/))
            case(252_i2)   !ST_int
                readStataVar(1:stataFile%nobs) = transfer(stataFile%theseData(varStart:(varStart+1),:),(/0_i2/))
            case(253_i2)   !ST_long
                readStataVar(1:stataFile%nobs) = transfer(stataFile%theseData(varStart:(varStart+3),:),(/0_i4/))
            case(254_i2)   !ST_float
                readStataVar(1:stataFile%nobs) = transfer(stataFile%theseData(varStart:(varStart+3),:),(/0.0_sp/))
            case(255_i2)   !ST_double
                readStataVar(1:stataFile%nobs) = transfer(stataFile%theseData(varStart:(varStart+7),:),(/0.0_dp/))
            case default
                call statamodWarn('unknown data type; setting array to zero')
                readStataVar(1:stataFile%nobs) = 0
        end select
    !end if

else !disk access... can be slow but doesn't require as much memory

    !move file pointer to start of data

    filePos = stataFile%stataPosition%varpos(thisVarno)

    select case(stataFile%typlist(thisVarno))
        case(1_i2:244_i2) !str
            call statamodWarn('strings are not supported; setting array to zero')
            readStataVar = 0
        case(251_i2)   !ST_byte
            do i=1, min(sizeReadStata,stataFile%nobs)
                read(stataFile%unit,POS=filePos) tempST_byte
                readStataVar(i:i)  = tempST_byte
                filePos = filePos + stataFile%stataPosition%nvarOffset
            end do
        case(252_i2)   !ST_int
            do i=1, min(sizeReadStata,stataFile%nobs)
                read(stataFile%unit,POS=filePos) tempST_int
                readStataVar(i:i)  = tempST_int
                filePos = filePos + stataFile%stataPosition%nvarOffset
            end do
        case(253_i2)   !ST_long
            do i=1, min(sizeReadStata,stataFile%nobs)
                read(stataFile%unit,POS=filePos) tempST_long
                readStataVar(i:i)  = tempST_long
                filePos = filePos + stataFile%stataPosition%nvarOffset
            end do
        case(254_i2)   !ST_float
            do i=1, min(sizeReadStata,stataFile%nobs)
                read(stataFile%unit,POS=filePos) tempST_float
                readStataVar(i:i)  = tempST_float
                filePos = filePos + stataFile%stataPosition%nvarOffset
            end do
        case(255_i2)   !ST_double
            do i=1, min(sizeReadStata,stataFile%nobs)
                read(stataFile%unit,POS=filePos) tempST_double
                readStataVar(i:i)  = tempST_double
                filePos = filePos + stataFile%stataPosition%nvarOffset
            end do
        case default
            call statamodWarn('unknown data type; setting array to zero')
            readStataVar = 0
    end select

end if

if (stataFile%doSample) then
    readStataVar2 = readStataVar
    do i = 1, stataFile%nobs
        readStataVar(i) = readStataVar2(stataFile%sampleWeight(i))
    end do
end if

