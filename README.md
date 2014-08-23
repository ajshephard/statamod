
STATAMOD
========

STATAMOD is a Fortran module that provides read/write access for Stata datasets 
from within Fortran. It is written by Andrew Shephard (<asheph@econ.upenn.edu>)


Stata read support
------------------

The reading of variables from a Stata .dta file is supported through the following 
subroutines and functions.

```subroutine openStata(filename [,cache])```

This opens the Stata .dta file *filename*. A dataset must be open before any other
(reading) STATAMOD commands can be used. By default, the entire dataset is loaded
into memory (`cache=.true.`). If you specify `openStata(filename, cache=.false.)` 
then it will not do this. Instead, it will just load the variables as requested 
from disk. As variables are stored non-sequentially in the Stata .dta format, this 
can be significantly slower than memory access.

```subroutine descStata()```

This describes the open dataset, printing output to the standard output unit. It is
much like Stata's own describe command.

`integer function nobsStata()`

Returns an integer equal to the number of observations in the open dataset.

`integer function nvarStata()`

Returns an integer equal to the number of variables in the open dataset.

`logical function existStata(varname)`

Returns *.true.* or *.false.* depending on whether the variable exists in the open
dataset.

`subroutine readStata(readStataVar,varname)`

Copies data from the Stata variable varname in the open dataset to *readStataVar*.
If *nobsStata* exceeds the dimension of *readStataVar* an error is raised.
Otherwise the data is copied to first *nobsStata* elements of *readStataVar*.

`subroutine closeOpenStata()`

Closes the open Stata dataset. Call this procedure after all the relevant data
has been read. This will close the file connected to the specified unit or free
memory used to temporarily store the dataset depending upon whether disk or
memory (default) mode was specified in the `openStata()` subroutine.

Stata read support example
----
An example of typical STATAMOD usage is provided below.

```
use statamod

implicit none

integer, parameter :: dp = selected_real_kind(15,100)
integer :: nObs

real(dp),  dimension(:),   allocatable :: hhincome
real(dp),  dimension(:,:), allocatable :: hours

call openStata('statafile.dta')
call descStata()

nObs = nobsStata()
allocate(hhincome(nObs))
allocate(hours(nObs,2))

call readStata(hhincome,'income')
call readStata(hours(:,1),'hours_male')
call readStata(hours(:,2),'hours_female')

call closeOpenStata()
```

Stata write support
-------------------

There are a small number of subroutines that provide Stata .dta write functionality.

`subroutine saveStata(fileName, obs [,label])`
            
This opens the Stata .dta file filename for saving. You must specify the 
number of observations *obs* to save (saving is on a variable by variable basis), 
and you can optionaly specify a label for the dataset. It will save the dataset 
in the Stata 7 SE file format.

`subroutine writeStata(thisVar, thisName [,thisLabel])`

Saves the Fortran variable *thisVar* with the Stata variable name *thisName*. You can 
optionally label the variable with *thisLabel*. You should call `writeStata` for 
every variable you wish to save.

`subroutine closeSaveStata()`

Closes the open saveStata dataset. The file will only be written when this
routine is called.


Further information
-------------------

If you receive a stack overflow error when using the subroutine readStata()
increase the stack size, or specify cache=.false. in the calling options.
This may happen with large data sets.

STATAMOD should be able to read Stata .dta files for versions 5 to 11. I have not
been able to test data created with Stata version 5 or 6 and neither do I have a
detailed description of the file format implemented in these versions. Please
contact me if such datasets produces any unexpected results

Software license
----------------

STATAMOD is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

STATAMOD is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with STATAMOD.  If not, see <http://www.gnu.org/licenses/>.
