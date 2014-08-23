
!    This file is part of STATAMOD. Copyright (c) 2006-2010 Andrew Shephard
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

allocate(stataFile%saveCurr%next, STAT=ios); 
                            
if (ios /= 0) then
    call statamodError('error allocating memory')
    !close(stataFile%saveUnit)
end if
                
stataFile%saveCurr => stataFile%saveCurr%next;
nullify(stataFile%saveCurr%next)
            
! keep track of variables
stataFile%saveNVar = stataFile%saveNVar + 1
