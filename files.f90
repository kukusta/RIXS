  SUBROUTINE ADD_FILE(fname,unit,alias)
    use m_files
    implicit none
    character(filename_length) :: fname
    character(aliasname_length), optional :: alias
    integer unit
! local vars
    character(256) msg
    integer iok
    type(file), pointer :: this
    type(file), target, allocatable :: newnode

    if(.not. allocated(headfile))then
      allocate(headfile,stat=iok,errmsg=msg)
      if(iok /= 0)call print_allocation_error('headfile',iok,msg)
      tail=>headfile
    else
      allocate(newnode,stat=iok,errmsg=msg)
      if(iok /= 0)call print_allocation_error('newnode',iok,msg)
      this=>headfile
      do
        if(.not. associated(this%next))then
          this%next=>newnode
          newnode%prev=>this
          this=>newnode
          exit
        endif
        this=>this%next
      enddo
    endif
    write(*,*)'this%alias=',this%alias
    this%fname(:)=fname(:)
    this%unit=unit
    if(present(alias))this%alias(:)=alias(:)

  END SUBROUTINE ADD_FILE

  INTEGER FUNCTION GET_UNIT(alias)
    use m_files
    implicit none
    character(aliasname_length) alias
    type(file), pointer :: this

    this=>headfile
    do
      this=>this%next
      if(.not. associated(this))then
        write(*,*)'  cannot find file with alias '//trim(alias)
        get_unit=0
        return
      endif
    enddo
  END FUNCTION GET_UNIT

  SUBROUTINE PRINT_FILES
    use m_files
    implicit none
    type(file), pointer :: this

    this=>headfile
    do while(associated(this))
      write(*,*)
      write(*,*)'unit= ',this%unit,' alias= ',trim(this%alias)&
           ,' filename= ',trim(this%fname)
      if(.not. associated(this%next))return
      this=>this%next
    enddo
  END SUBROUTINE PRINT_FILES

  SUBROUTINE REMOVE_FILES
    use m_files
  END SUBROUTINE REMOVE_FILES
