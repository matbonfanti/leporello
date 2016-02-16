!***************************************************************************************
!*                              MODULE PrintTools
!***************************************************************************************
!
!>  \brief         Plot the data in VTK XML format
!>  \details       The following operations are performed with the subroutine
!>                 of this module:
!>      \arg  INITIALIZE A COLLECTION OF VTK FILES ( VTK_NewCollection ) \n
!>            In this initialization an array is allocated to store the names of
!>            the files composing the collection and the name of the collection is
!>            defined. This step is not mandatory. In case a non-initialized collection
!>            is used, it will be assumed that the user want to print a single file.
!>      \arg  SET A NEW FILE WITH A RECTILINEAR OR STRUCTURED GRID \n
!>            Given the name of the file and the mesh data, a new VTK file is started.
!>            The file is open and the grid is written in the appropriate format: 
!>            rectilinear or structured grid. 
!>            The name is optional: if it is not given, it will be used the name of the
!>            collection plus the ordinal number of the file in the collection. 
!>      \arg  WRITE SCALAR FIELD \n
!>            Write a scalar field in a file that already has a defined mesh. By default,
!>            the file is closed. If required, the file can be let open to add more than one
!>            scalar field to the same mesh. 
!>      \arg  WRITE COLLECTION FILE \n
!>            when at least one file has been closed and if no VTK file is currently in
!>            output, the collection file can be written. This file define the VTK files 
!>            in the collection as time snapshots.
!>      \arg  WRITE TRAJECTORY \n
!>            instead of rectilinear/structured grid + data, collection slots can be used
!>            to write trajectory files. A trajectory file is a VTP (polydata) file 
!>            specifying the positions in time and line connection between the positions 
!
!***************************************************************************************
!
!>  \author        Matteo Bonfanti
!>  \version       2.0
!>  \date          May 2010
!>
!***************************************************************************************
!
!>  \par Updates
!>  \arg The function LookForFreeUnit has been moved to the module MyConsts
!>  \arg Major revision: completely different organization of snapshot and
!>                       collections (12 July 2012)
!>  \arg Implemented the possibility to define structured grid (13 July 2012) 
!>  \arg Implemented and tested the possibility to write trajectory VTP files (13 July 2012) 
!
!***************************************************************************************
!
!>  \todo Add the possibility to print data in binary form through preprocessing directives
!>  \todo Implement the possibility of including vector field in the file
!
!***************************************************************************************
MODULE PrintTools
#include "preprocessoptions.cpp"

   IMPLICIT NONE

   PRIVATE

   ! data type for VTK collection and files
   PUBLIC  :: VTKInfo
   ! initialize and print a collection of VTK files
   PUBLIC  :: VTK_NewCollection, VTK_PrintCollection
   ! define a new file with a given mesh (structured grid, rectilinear grid )
   PUBLIC  :: VTK_NewStructuredSnapshot, VTK_NewRectilinearSnapshot
   ! add a field to the vtk file with mesh
   PUBLIC  :: VTK_AddScalarField
   ! write a 1D or 2D or 3D trajectory to a polydata VTK file
   PUBLIC  :: VTK_WriteTrajectory

   ! Data formats
   ! strings with increasing and decreasing indentation
   CHARACTER(10), PARAMETER, PRIVATE :: F1 = "(100A)    "
   CHARACTER(10), PARAMETER, PRIVATE :: F2 = "(2X,98A)  "
   CHARACTER(10), PARAMETER, PRIVATE :: F3 = "(4X,96A)  "
   CHARACTER(10), PARAMETER, PRIVATE :: F4 = "(6X,94A)  "
   CHARACTER(10), PARAMETER, PRIVATE :: F5 = "(8X,92A)  "
   CHARACTER(10), PARAMETER, PRIVATE :: F6 = "(10X,90A) "
   ! numeric data on six columns
   CHARACTER(20), PARAMETER, PRIVATE :: FN = "(1(10X,6(E15.8,X)))"
   CHARACTER(20), PARAMETER, PRIVATE :: FI = "(1(10X,18(I5,X)))"

   ! status of the collection
   INTEGER, PARAMETER, PRIVATE :: NON_INITIALIZED = 0
   INTEGER, PARAMETER, PRIVATE :: INITIALIZED     = 1
   INTEGER, PARAMETER, PRIVATE :: MESH_WRITTEN    = 2
   INTEGER, PARAMETER, PRIVATE :: DATA_WRITTEN    = 3

   ! kind of mesh written for the current file
   INTEGER, PARAMETER, PRIVATE :: NO_MESH         = 0
   INTEGER, PARAMETER, PRIVATE :: RECTILINEAR     = 1
   INTEGER, PARAMETER, PRIVATE :: STRUCTURED      = 2
   
   !> VTK files and collection data type
   TYPE :: VTKInfo
      PRIVATE
      INTEGER    :: Status = NON_INITIALIZED                     !< Actual status of the VTK type
      CHARACTER(len=30)    :: CollectionName                     !< The name of the collection
      CHARACTER(len=30), DIMENSION(:), ALLOCATABLE :: FileNames  !< The name of the files
      INTEGER    ::  Counter = 0                                 !< The VTK file that is being written
      INTEGER    ::  Unit                                        !< The unit to which it is being written
      INTEGER    ::  ni, nj, nk                                  !< The mesh size of the current VTK file
      INTEGER    ::  MeshKind = NO_MESH                          !< kind of mesh written for the current file
   END TYPE VTKInfo

   
!********************************************************************************************************
   CONTAINS
!********************************************************************************************************



!*******************************************************************************
!          VTK_NewCollection
!*******************************************************************************
!> Initialize a VTK file collection.
!> In detail, the subroutine allocate an array where to store the names of 
!> the individual VTK files, and set the name of the collection
!>
!> @param VTKType            VTK collection data type
!> @param CollectionSize     Number of maximum file in the collection
!> @param CollectionName     Name of the collection
!*******************************************************************************
   SUBROUTINE VTK_NewCollection ( VTKType, CollectionSize, CollectionName )
      IMPLICIT NONE
      TYPE(VTKInfo), INTENT(inout)      :: VTKType
      INTEGER, INTENT(in)               :: CollectionSize
      CHARACTER(len=*), INTENT(in)      :: CollectionName

      ! Check if the collection has been already initialized
      CALL WARN( VTKType%Status /= NON_INITIALIZED, "VTK_NewCollection: VTKType already initialized" &
                  // NewLine // " Previous data will be overwritten " )

      ! Store collection name
      VTKType%CollectionName = TRIM(ADJUSTL( CollectionName ))

      ! Allocate filename array
      ALLOCATE( VTKType%FileNames( CollectionSize ) )

      ! Set actual status 
      VTKType%Status = INITIALIZED
      
   END SUBROUTINE VTK_NewCollection


!*******************************************************************************
!          VTK_NewRectilinearSnapshot
!*******************************************************************************
!> Initialize a new VTR file and set the rectilinear mesh.
!> If the collection was not initialize, do it with for just 1 file.
!> The name of the file is optional. If not given, the name of the collection
!> is used. Since many input variables are optional, it may be necessary to 
!> explicitly identify them calling the subroutine:
!> CALL VTK_NewRectilinearSnapshot ( VTK, X=.., Y=.., FileName=".." )
!>
!> @param VTKType            VTK collection data type
!> @param X                  Real array with the X of the grid points
!> @param Y                  Real array with the Y of the grid points (optional)
!> @param Z                  Real array with the Z of the grid points (optional)
!> @param FileName           Name of the file (optional)
!*******************************************************************************
   SUBROUTINE VTK_NewRectilinearSnapshot ( VTKType, X, Y, Z, FileName )
      IMPLICIT NONE
      TYPE(VTKInfo), INTENT(inout)              :: VTKType
      REAL, DIMENSION(:), INTENT(in)            :: X
      REAL, DIMENSION(:), INTENT(in), OPTIONAL  :: Y, Z
      CHARACTER(len=*), OPTIONAL, INTENT(in)    :: FileName
      INTEGER                                   :: ierr, u
      CHARACTER(len=40)                         :: MeshSize
      CHARACTER(len=10)                         :: Dim1, Dim2, Dim3
      
      ! Check and in case initialize the collection
      IF (VTKType%Status == NON_INITIALIZED)  CALL VTK_NewCollection( VTKType, 1, "Unknown" )
      
      ! Give error if the previous file was not completed
      CALL ERROR( VTKType%Status == MESH_WRITTEN .OR. VTKType%Status == DATA_WRITTEN, &
                    "VTK_NewRectilinearSnapshot: trying to open a new snapshot without closing the previous one" )

      ! Find free I/O unit
      VTKType%Unit      = LookForFreeUnit()      
      ! Increase counter
      VTKType%Counter   = VTKType%Counter + 1                    

      ! Check if the number of available file in the collection is over
      CALL ERROR( VTKType%Counter > size(VTKType%FileNames), &
               "VTK_NewRectilinearSnapshot: no more file available in the VTK collection" )
      
      ! Define snapshot name
      IF ( PRESENT( FileName ) ) THEN                             ! in case the namefile has been given
            VTKType%FileNames(VTKType%Counter) = TRIM(ADJUSTL( FileName )) // ".vtr"
      ELSE                                                        ! otherwise the collection name is used
            VTKType%FileNames(VTKType%Counter) = TRIM(VTKType%CollectionName) // "_" // &
                                                 NumberToString( VTKType%Counter ) // ".vtr"
      END IF

      ! Define mesh sizes
      VTKType%ni=size(X)
      IF ( Present(Y) ) THEN
               VTKType%nj = size(Y)
      ELSE
               VTKType%nj = 1
      ENDIF
      IF ( Present(Z) ) THEN
               VTKType%nk = size(Z)
      ELSE
               VTKType%nk = 1
      ENDIF
      
      ! Build string with numbers of grid size
      WRITE( Dim1 , '(I6)') VTKType%ni
      WRITE( Dim2 , '(I6)') VTKType%nj
      WRITE( Dim3 , '(I6)') VTKType%nk
      WRITE( MeshSize, '("1 ",A," 1 ",A," 1 ",A)' ) TRIM(ADJUSTL(Dim1)), TRIM(ADJUSTL(Dim2)), TRIM(ADJUSTL(Dim3))
      
      ! Open file unit
      u = VTKType%Unit
      OPEN( UNIT=u, FILE=TRIM(VTKType%FileNames(VTKType%Counter)), FORM="formatted", STATUS="replace", ACTION="write", IOSTAT=ierr )
      CALL ERROR( (ierr/=0), "VTK_NewRectilinearSnapshot: Error creating file" // TRIM(VTKType%FileNames(VTKType%Counter)) )

      WRITE(u,F1)   '<VTKFile type="RectilinearGrid" format="ascii">'
      WRITE(u,F2)   '<RectilinearGrid WholeExtent=" ' // TRIM(MeshSize) // ' ">'
      WRITE(u,F3)   '<Piece Extent=" ' // TRIM(MeshSize) // ' ">'
      WRITE(u,F4)   '<Coordinates>'
      WRITE(u,F5)   '<DataArray type="Float64" Name="X_COORDINATES" NumberOfComponents="1" format="ascii">'
      WRITE(u,FN)        X( 1 : VTKType%ni )
      WRITE(u,F5)   '</DataArray>'
      WRITE(u,F5)   '<DataArray type="Float64" Name="Y_COORDINATES" NumberOfComponents="1" format="ascii">'
      IF (Present(Y)) THEN
         WRITE(u,FN)     Y( 1 : VTKType%nj )
      ELSE
         WRITE(u,FN)     0.0000
      ENDIF
      WRITE(u,F5)   '</DataArray>'
      WRITE(u,F5)   '<DataArray type="Float64" Name="Z_COORDINATES" NumberOfComponents="1" format="ascii">'
      IF (Present(Z)) THEN
         WRITE(u,FN)     Z( 1 : VTKType%nk )
      ELSE
         WRITE(u,FN)     0.0000
      ENDIF
      WRITE(u,F5)   '</DataArray>'
      WRITE(u,F4)   '</Coordinates>'
      WRITE(u,F4)   '<PointData>'      
      
      ! Set the status to mesh_written
      VTKType%Status = MESH_WRITTEN

      ! Define the type of mesh written
      VTKType%MeshKind = RECTILINEAR
      
   END SUBROUTINE VTK_NewRectilinearSnapshot

!*******************************************************************************
!          VTK_NewStructuredSnapshot
!*******************************************************************************
!> Initialize a new VTS file and set the structured grid.
!> If the collection was not initialize, do it with for just 1 file.
!> The name of the file is optional. If not given, the name of the collection
!> is used. Since many input variables are optional, it may be necessary to 
!> explicitly identify them calling the subroutine:
!> CALL VTK_NewRectilinearSnapshot ( VTK, X=.., Y=.., FileName=".." )
!>
!> @param VTKType            VTK collection data type
!> @param Coords             Array 3xN with the coordinates of the N grid points
!> @param N1                 First dimension of the grid
!> @param N2                 Second dimension of the grid (optional, if not given 1 is assumed)
!> @param N3                 Third dimension of the grid (optional, if not given 1 is assumed)
!> @param FileName           Name of the file (optional)
!*******************************************************************************
   SUBROUTINE VTK_NewStructuredSnapshot ( VTKType, Coords, N1, N2, N3, FileName )
      IMPLICIT NONE
      TYPE(VTKInfo), INTENT(inout)              :: VTKType
      REAL, DIMENSION(:,:), INTENT(in)          :: Coords
      INTEGER, INTENT(in)                       :: N1
      INTEGER, INTENT(in), OPTIONAL             :: N2, N3
      CHARACTER(len=*), OPTIONAL, INTENT(in)    :: FileName
      INTEGER                                   :: ierr, u, i
      CHARACTER(len=40)                         :: MeshSize
      CHARACTER(len=10)                         :: Dim1, Dim2, Dim3
      
      ! Check and in case initialize the collection
      IF (VTKType%Status == NON_INITIALIZED)  CALL VTK_NewCollection( VTKType, 1, "Unknown" )
      
      ! Give error if the previous file was not completed
      CALL ERROR( VTKType%Status == MESH_WRITTEN .OR. VTKType%Status == DATA_WRITTEN, &
                    "VTK_NewStructuredSnapshot: trying to open a new snapshot without closing the previous one" )

      ! Find free I/O unit
      VTKType%Unit      = LookForFreeUnit()      
      ! Increase counter
      VTKType%Counter   = VTKType%Counter + 1                    

      ! Check if the number of available file in the collection is over
      CALL ERROR( VTKType%Counter > size(VTKType%FileNames), &
               "VTK_NewStructuredSnapshot: no more file available in the VTK collection" )
      
      ! Define snapshot name
      IF ( PRESENT( FileName ) ) THEN                             ! in case the namefile has been given
            VTKType%FileNames(VTKType%Counter) = TRIM(ADJUSTL( FileName )) // ".vts"
      ELSE                                                        ! otherwise the collection name is used
            VTKType%FileNames(VTKType%Counter) = TRIM(VTKType%CollectionName) // "_" // &
                                                 NumberToString( VTKType%Counter ) // ".vts"
      END IF

      ! Define mesh sizes
      VTKType%ni= N1
      IF ( Present(N2) ) THEN
               VTKType%nj = N2
      ELSE
               VTKType%nj = 1
      ENDIF
      IF ( Present(N3) ) THEN
               VTKType%nk = N3
      ELSE
               VTKType%nk = 1
      ENDIF

      ! check the size of the structured grid
      CALL ERROR ( size(Coords(:,1)) /= 3 ,  "VTK_NewStructuredSnapshot : wrong coordinates array dimension (1)." )
      CALL ERROR ( size(Coords(1,:)) /= VTKType%ni * VTKType%nj * VTKType%nk ,    &
             "VTK_NewStructuredSnapshot : wrong coordinates array dimension (2)." )
      
      ! Build string with numbers of grid size
      WRITE( Dim1 , '(I6)') VTKType%ni
      WRITE( Dim2 , '(I6)') VTKType%nj
      WRITE( Dim3 , '(I6)') VTKType%nk
      WRITE( MeshSize, '("1 ",A," 1 ",A," 1 ",A)' ) TRIM(ADJUSTL(Dim1)), TRIM(ADJUSTL(Dim2)), TRIM(ADJUSTL(Dim3))
      
      ! Open file unit
      u = VTKType%Unit
      OPEN( UNIT=u, FILE=TRIM(VTKType%FileNames(VTKType%Counter)), FORM="formatted", STATUS="replace", ACTION="write", IOSTAT=ierr )
      CALL ERROR( (ierr/=0), "VTK_NewStructuredSnapshot: Error creating file" // TRIM(VTKType%FileNames(VTKType%Counter)) )

      WRITE(u,F1)   '<VTKFile type="StructuredGrid" format="ascii">'
      WRITE(u,F2)   '<StructuredGrid WholeExtent=" ' // TRIM(MeshSize) // ' ">'
      WRITE(u,F3)   '<Piece Extent=" ' // TRIM(MeshSize) // ' ">'
      WRITE(u,F4)   '<Points>'
      WRITE(u,F5)   '<DataArray type="Float64" Name="GRID_COORDS" NumberOfComponents="3" format="ascii">'
      WRITE(u,FN)    ( Coords(1:3, i), i = 1 , VTKType%ni*VTKType%nj*VTKType%nk )
      WRITE(u,F5)   '</DataArray>'
      WRITE(u,F4)   '</Points>'
      WRITE(u,F4)   '<PointData>'      
      
      ! Set the status to mesh_written
      VTKType%Status = MESH_WRITTEN

      ! Define the type of mesh written
      VTKType%MeshKind = STRUCTURED
      
   END SUBROUTINE VTK_NewStructuredSnapshot


!*******************************************************************************
!          VTK_AddScalarField
!*******************************************************************************
!> Add a scalar field to an open VTK file.
!> By default, the subroutine closes the VTK file. If otherwise specified (by
!> means of the variable LetFileOpen) the closing strings are not added and the
!> file is not closed.
!>
!> @param VTKType            VTK collection data type
!> @param Name               Name of the field
!> @param Field              Real array with the field values
!> @param LetFileOpen        Do not close the file (optional argument) 
!*******************************************************************************
   SUBROUTINE VTK_AddScalarField (VTKType, Name, Field, LetFileOpen )
      IMPLICIT NONE
      TYPE(VTKInfo), INTENT(inout)              :: VTKType
      CHARACTER(len=*), INTENT(in)              :: Name
      REAL, INTENT(in), DIMENSION(:)            :: Field
      LOGICAL, INTENT(in), OPTIONAL             :: LetFileOpen
      INTEGER                                   :: u
      LOGICAL                                   :: CloseFile

      ! check the status of the VTK data type
      CALL ERROR( (VTKType%Status /= MESH_WRITTEN) .AND. (VTKType%Status /= DATA_WRITTEN), &
            "VTK_AddScalarField : trying to write scalar field to a non initialized file." )

      ! check the size of the scalar field
      CALL WARN( (size(Field)/=VTKType%ni*VTKType%nj*VTKType%nk),    &
            "VTK_AddScalarField : Incompatible FIELD and MESH sizes." )

      ! dummy unit variable
      u = VTKType%Unit
      
      ! write scalar field
      WRITE(u,F5) '<DataArray type="Float64" Name="'//TRIM(ADJUSTL(Name))//'" NumberOfComponents="1" format="ascii">'
      WRITE(u,FN) Field(1:VTKType%ni*VTKType%nj*VTKType%nk)
      WRITE(u,F5) '</DataArray>'
      
      ! set the status to data written
      VTKType%Status = DATA_WRITTEN
     
      ! decide whether to close the file
      IF ( PRESENT(LetFileOpen) ) THEN
           CloseFile = .NOT. LetFileOpen
      ELSE
           CloseFile = .TRUE.
      END IF

      ! close the file
      IF ( CloseFile ) THEN
         ! add last lines
         WRITE(u,F4) '</PointData>'
         WRITE(u,F3) '</Piece>'
         IF ( VTKType%MeshKind == RECTILINEAR ) THEN
               WRITE(u,F2) '</RectilinearGrid>'
         ELSE IF ( VTKType%MeshKind == STRUCTURED ) THEN
               WRITE(u,F2) '</StructuredGrid>'
         ENDIF
         WRITE(u,F1) '</VTKFile>'
         ! close the file
         CLOSE( u )
         ! set the status to initialized: new file can be started 
         VTKType%Status = INITIALIZED
      END IF
      
   END SUBROUTINE VTK_AddScalarField

!*******************************************************************************
!          VTK_WriteTrajectory
!*******************************************************************************
!> Print a trajectory in VTP format. 
!> A 1D, 2D or 3D trajectory can be provided, with an array NDim times NTstep
!> where NDim is the number of dimension, and NTstep is the arbitrary number 
!> of time steps. The trajectory is printed in a "PolyData" type file at once.
!> The file is closed and nothing else can be added.
!>
!> @param VTKType            VTK collection data type
!> @param TrajXYZ            Array with the trajectory to be printed
!> @param FileName           Name of the file (optional)
!*******************************************************************************
  SUBROUTINE VTK_WriteTrajectory ( VTKType, TrajXYZ, FileName )
    IMPLICIT NONE
      TYPE(VTKInfo), INTENT(inout)              :: VTKType
      REAL, DIMENSION(:,:), INTENT(IN)          :: TrajXYZ
      CHARACTER(len=*), OPTIONAL, INTENT(in)    :: FileName
      INTEGER                                   :: u, NPts, i
      CHARACTER(10)                             :: No1, No2

      ! Check and in case initialize the collection
      IF (VTKType%Status == NON_INITIALIZED)  CALL VTK_NewCollection( VTKType, 1, "Unknown" )
      
      ! Give error if the previous file was not completed
      CALL ERROR( VTKType%Status == MESH_WRITTEN .OR. VTKType%Status == DATA_WRITTEN, &
                    "VTK_WriteTrajectory: trying to open a new snapshot without closing the previous one" )

      ! Find free I/O unit
      VTKType%Unit      = LookForFreeUnit()      
      ! Increase counter
      VTKType%Counter   = VTKType%Counter + 1    

      ! Check if the number of available file in the collection is over
      CALL ERROR( VTKType%Counter > size(VTKType%FileNames), &
               "VTK_WriteTrajectory: no more file available in the VTK collection" )

      ! Define snapshot name
      IF ( PRESENT( FileName ) ) THEN                             ! in case the namefile has been given
            VTKType%FileNames(VTKType%Counter) = TRIM(ADJUSTL( FileName )) // ".vtp"
      ELSE                                                        ! otherwise the collection name is used
            VTKType%FileNames(VTKType%Counter) = TRIM(VTKType%CollectionName) // "_" // &
                                                 NumberToString( VTKType%Counter ) // ".vtp"
      END IF

      ! check the size of the structured grid
      CALL ERROR ( size(TrajXYZ(:,1)) > 3 ,  "VTK_WriteTrajectory : wrong coordinates array dimension (1)." )

      ! Store the other dimension
      NPts = size(TrajXYZ(1,:))

      ! Build string with numbers of grid size
      WRITE( No1 , '(I10)') NPts
      WRITE( No2 , '(I10)') NPts-1

      ! Open file unit
      u = VTKType%Unit
      OPEN( UNIT=u, FILE=TRIM(VTKType%FileNames(VTKType%Counter)), FORM="formatted", STATUS="replace", ACTION="write", IOSTAT=i )
      CALL ERROR( (i/=0), "VTK_WriteTrajectory: Error creating file" // TRIM(VTKType%FileNames(VTKType%Counter)) )

      WRITE(u,F1)   '<VTKFile type="PolyData" format="ascii">'
      WRITE(u,F2)   '<PolyData>'
      WRITE(u,F3)   '<Piece NumberOfPoints="'//TRIM(ADJUSTL(No1))// '" NumberOfVerts="0" NumberOfLines="' &
                         // TRIM(ADJUSTL(No2)) // '" NumberOfStrips="0" NumberOfPolys="0">'

      ! Write trajectory points
      WRITE(u,F4)   '<Points>'
      WRITE(u,F5)   '<DataArray type="Float64" NumberOfComponents="3" format="ascii">'
      SELECT CASE ( size(TrajXYZ(:,1)) )
         CASE ( 1 )
            WRITE(u,FN)    ( (/ TrajXYZ(1,i), 0.0, 0.0 /), i = 1, NPts )
         CASE ( 2 )
            WRITE(u,FN)    ( (/ TrajXYZ(1,i), TrajXYZ(2,i), 0.0 /), i = 1, NPts )
         CASE ( 3 )
            WRITE(u,FN)    ( TrajXYZ(1:3,i), i = 1, NPts )
      END SELECT
      WRITE(u,F5)   '</DataArray>'
      WRITE(u,F4)   '</Points>'

      ! Connect trajectory points with lines
      WRITE(u,F4)   '<Lines>'
      WRITE(u,F5) '<DataArray type="Float64" Name="connectivity" NumberOfComponents="1" format="ascii">'
      WRITE(u,FI)  ( (/ i-1, i /), i = 1,NPts-1 )
      WRITE(u,F5) '</DataArray>'
      WRITE(u,F5) '<DataArray type="Float64" Name="offsets" NumberOfComponents="1" format="ascii">'
      WRITE(u,FI)  ( 2*i, i = 1,NPts-1 )
      WRITE(u,F5) '</DataArray>'
      WRITE(u,F4)   '</Lines>'      

      WRITE(u,F3) '</Piece>'
      WRITE(u,F2) '</PolyData>'
      WRITE(u,F1) '</VTKFile>'

      ! close the file
      CLOSE( u )
      ! set the status to initialized: new file can be started 
      VTKType%Status = INITIALIZED

  END SUBROUTINE VTK_WriteTrajectory


!*******************************************************************************
!          VTK_PrintCollection
!*******************************************************************************
!> Print the collection XML file. 
!>
!> @param VTKType            VTK collection data type
!*******************************************************************************
   SUBROUTINE VTK_PrintCollection ( VTKType )
      IMPLICIT NONE
      TYPE(VTKInfo), INTENT(inout) :: VTKType
      INTEGER                      :: u, ierr, Shot
      CHARACTER(10)                :: SnapNumber
      
      ! to print the collection, at least one file has to be written and
      ! the last file has to be closed
      CALL ERROR( VTKType%Counter < 1, "VTK_PrintCollection : not even one file in the collection ") 
      CALL ERROR( VTKType%Status /= INITIALIZED, "VTK_PrintCollection : Last file is still to be closed ") 

      ! find a free unit to write the collection file
      u = LookForFreeUnit()

      ! open unit
      OPEN( UNIT=u, FILE=TRIM(VTKType%CollectionName)//".pvd",  &
            FORM="FORMATTED", STATUS="replace", ACTION="write", IOSTAT=ierr)
      CALL ERROR( (ierr/=0), "Error creating file"//TRIM(VTKType%CollectionName)//".pvd" )
      
      ! write first lines 
      WRITE(u,F1)    '<?xml version="1.0"?>'
      WRITE(u,F1)    '<VTKFile type="Collection" format="ascii">'
      WRITE(u,F2)    '<Collection>'

      ! write file list with increasing time
      DO Shot = 1, VTKType%Counter
            WRITE( SnapNumber,"(I10)") (Shot-1)
            WRITE(u,F3)  '<DataSet timestep="' // TRIM(ADJUSTL(SnapNumber)) // '" part="0" file="' &
                                                                 // TRIM(VTKType%FileNames(Shot)) //'"/>'
      END DO

      ! write last lines
      WRITE(u, F2)   '</Collection>'
      WRITE(u, F1)   '</VTKFile>'
      
      ! close the file
      CLOSE( u )

   END SUBROUTINE VTK_PrintCollection
   
END MODULE PrintTools
