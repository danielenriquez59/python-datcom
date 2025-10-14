

THE USAF DIGITAL DATCOM                                    /datcom/readme.txt

The files for this program are in the directory datcom 
  readme.txt      general description
  namelist.pdf    descriptions of the namelist data items
  datcom.f        The complete source code
  exwin.zip       The 11 sample cases from AFFDL-TR-79-3032 along with the
                    expected output for each. Input files have the Window/DOS
                    end-of-line characters 
  exlinux.zip     Same as exwin.zip, but with Unix end-of line

The reference documents for Digital Datcom may be accessed
from the web page https://www.pdas.com/datcomrefs.html. The description
of the input to Digital Datcom is in AFFDL-TR-79-3032.
 
This program asks for the name of the input file. This must be 
a file written in the input format to Digital Datcom. After reading 
the input data, the program produces a file called datcom.out that is
formatted for printing. The printed output is designed for printers with
132 columns, so a modern user must set the printer to landscape mode and
8-point type to keep the print on the page.

To compile this program for your machine, use the command
   gfortran  datcom.f  -o datcom.exe
Linux and Macintosh users may prefer to omit the .exe on the file name.
You will get hundreds, even thousands of warning messages when compiling
as this program does not conform to modern usage. However, they are all
warnings; there should be no errors.


Years ago, the USAF commissioned a project devoted to the idea of
documenting methods for estimating stability and control
characteristics of airplanes and missiles. The result was a big four-
notebook document. That document is the basis for this program. 
The entire document is available for download as a PDF from the 
references page at https://www.pdas.com/datcomrefs.html.
Before you decide to print the entire document yourself, you should
work out your cost in paper and cartridges and time to print 3134 pages.
There is a firm called IHS (previously, Global Engineering Documents) that 
offers two-sided copies in loose-leaf binders. 
As of 18 July 2013, their price was $195.
You can send e-mail to global@ihs.com or visit them at 
https://global.ihs.com   and purchase online.

A natural development of this project with its printed charts and tables
was a computer program that would do the table lookup and interpolation.
The result was the program known as digital datcom.

The manual for Digital Datcom is quite different from the
Datcom theory document. It refers to the usage of the computer program.


The DATCOM cases use $ as the beginning and ending namelist characters
instead of the modern usage. If you look at the code carefully, you 
will see that the native namelist of Fortran is not being used, but 
that an entire namelist emulator has been written for Datcom. This 
was probably a smart move back in the 60s and 70s before namelist 
was standardized.

There are 11 sample cases (both input and output) supplied so you can
be sure your system is operating correctly. However, if you use the
files ex1.dat, etc. on either Linux or Macintosh, you will probably 
get an error message indicating that there is an error in the input
file. This happens because the Windows end-of-line characters are not
the same, and there appears to be an extraneous illegal character.
There is a file called exlinux.zip that have the same
example cases with the proper end of line characters.
