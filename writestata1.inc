
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

character(*), intent(in)           :: thisName
character(*), intent(in), optional :: thisLabel

integer(i4)   :: i, ios
integer(i1)   :: tempInt1
character(33) :: tempStr33

if (.not. stataFile%saveInit) then
    call statamodError('a file has not been specified to write to')
end if

if (len(thisName)>32) then
    call statamodError('variable '//thisName//' exceeds 32 characters')
end if

write (TempStr33,*) stataFile%saveNObs

if (size(thisVar) /= stataFile%saveNObs) then
    call statamodError('variable ' // trim(thisName) // ' is the wrong dimension. '// &
        & 'Expecting an array of size ' // trim(adjustl(TempStr33)) // '*1')
end if

if (.not. stataFile%saveOnce) then

    ! first var to be written
    allocate(stataFile%saveCurr, STAT = ios)

    if (ios /= 0) then
        call statamodError('error allocating memory')
        !close(stataFile%saveUnit)
    end if

    if (associated(stataFile%saveCurr%next)) then
        nullify(stataFile%saveCurr%next)
    end if

    stataFile%saveHead => stataFile%saveCurr

    stataFile%saveOnce = .true.

end if

if (stataFile%saveNVar >= 32767) then
    call statamodError('attempt to exceed Stata maxvar limit (32767)')
    !close(stataFile%saveUnit)
end if

! var label
if (present(thisLabel)) then
    tempInt1 = len(thisLabel) + 1
    stataFile%saveCurr%label = thisLabel
    stataFile%saveCurr%label(tempInt1:81) = char(0)
else
    stataFile%saveCurr%label = char(0)
end if

! var name, first check whether it is stata legal
do i = 1, len(trim(thisName))
    select case(ichar(thisName(i:i)))
        case(48:57, 65:90, 95, 97:122) !0:9,a:z,_,A:Z
        case default
            call statamodError('variable name ' // thisName // ' contains an illegal character. ' // &
                & 'Only alpha-numeric characters and underscore are allowed')
    end select
end do !I

select case(ichar(thisName(1:1)))
    case(48:57) !0:9
        call statamodError('variable name ' // thisName // ' contains an illegal character. ' // &
            & 'Names can not begin with a numeric character')
    case default
end select

! now store in memory
tempInt1 = len(trim(thisName)) + 1
stataFile%saveCurr%name = trim(thisName)
stataFile%saveCurr%name(tempInt1:33) = char(0)

! and now, check unique
tempStr33 = stataFile%saveCurr%name
stataFile%saveCurr => stataFile%saveHead

do while(associated(stataFile%saveCurr%next))
    if (tempStr33 == stataFile%saveCurr%name) then
        call statamodError('a variable already exists with name ' // trim(tempStr33))
    end if
    stataFile%saveCurr => stataFile%saveCurr%next
end do
