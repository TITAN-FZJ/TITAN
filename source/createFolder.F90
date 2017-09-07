subroutine create_folder()
  use mod_parameters, only: itype
  use mod_SOC, only: SOCc
  implicit none


  character(len=500) :: base_command = ""
  character(len=500) :: folder = ""
  write(base_command, "('mkdir -p ./results/',a1,'SOC/')") SOCc

  !! Create selfconsistency folder
  write(folder, "(a,'selfconsistency')") trim(base_command)
  call execute_command_line(trim(folder))

  return
end subroutine create_folder
