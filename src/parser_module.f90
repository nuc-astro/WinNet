!> @file parser_module.f90
!!
!! The error file code for this file is ***W36***.
!!
!! @brief Module \ref parser_module, which enables parsing strings for analytic trajectory mode
!!

!> @brief Subroutines for equation parsing in case of an analytic trajectory
!!        or luminosity mode.
!!
!! The subroutine can calculate equations from a string, e.g., "5.1+3d-2*(2-1)*sin(pi)/x".
!! A variable name (in the upper example "x") can be given, which will be filled with a value.
!! Everything in this function is private with the exception of the function
!! "parse_string" which serves as a interface to parse the given string.
!! The return value of this function is a real (i.e., not a string).
!!
!!
!! @author Moritz Reichert
!! @date   20.12.20
!!
#include "macros.h"
module parser_module
   use parameter_class, only: max_fname_len
   use error_msg_class, only: raise_exception,num_to_str,int_to_str
      implicit none

      ! Parsing strings and lengths
      character(max_fname_len)              :: parsing_string    !< Current substring to parse
      character(max_fname_len)              :: complete_string   !< Complete input string

      ! Operator and number arrays that got parsed
      integer                                :: operator_count, number_count
      character(1),dimension(:),allocatable  :: operators  !< Array containing the operations of the equation as a character (e.g., ["*","+"])
      real(r_kind),dimension(:),allocatable  :: numbers    !< Array containing the numbers of the equation (e.g., [1, 2, 3])

      ! Variable to replace, the value is set by "parse_string"
      real(r_kind)                           :: variable_value    !< Input variable, set by \ref parse_string
      character(1),parameter                 :: variable_name='x' !< Name of the variable that is replaced by parse_string


     !
     ! Public and private fields and methods of the module
     ! Only the function parse_string and find_start_value
     ! should be public
     !
     public:: &
        parse_string,find_start_value
     private:: &
        parsing_string,complete_string,operator_count,number_count, &
        operators,numbers,variable_value, variable_name
     private:: &
        is_digit,is_operator,is_dot_operator,is_line_operator, &
        is_power_operator,is_separator,operation,eval_function, &
        make_pars_str_consistent,common_variables,store_digits_operator, &
        get_simple_result,evaluate

      contains


      !>
      !! @brief Function to decide whether a character is a number or not (i.e., 0-9).
      !!
      !! ### Example
      !!~~~~~~~~~~~~~~.f90
      !! b = is_digit("5")
      !! c = is_digit("a")
      !!~~~~~~~~~~~~~~
      !! b will be .True. and c will be .False.
      !!
      !!
      !! @author  Moritz Reichert
      !! @date    20.12.20
      function is_digit(str_in)
         implicit none
         character(1),intent(in)  :: str_in     !< Input character
         integer                  :: ascii_rep  !< ascii representation of string
         logical                  :: is_digit   !< True in case of a digit and False in case of another character

         ! Check numbers between 0 and 9 by their ascii representation
         ascii_rep = iachar(str_in(1:1))
         if (ascii_rep>= iachar("0") .and. ascii_rep<=iachar("9") ) then
            is_digit=.True.
         else
            is_digit = .False.
         endif
      end function is_digit


      !>
      !! @brief Function to decide whether a character is a mathematical operator
      !!        or not. Functions (as e.g., sin) are not defined here!
      !!
      !! ### Example
      !!~~~~~~~~~~~~~~.f90
      !! b = is_operator("+")
      !! c = is_operator("5")
      !!~~~~~~~~~~~~~~
      !! b will be .True. and c will be .False.
      !!
      !! @see is_line_operator, is_power_operator, is_separator, is_dot_operator
      !!
      !! @author  Moritz Reichert
      !! @date    20.12.20
      function is_operator(str_in)
         implicit none
         character(1),intent(in)  :: str_in       !< Input character
         logical                  :: is_operator  !< True in case of an operator (+,-,*,/,^) and False otherwise

         select case(str_in)
            case('+')
               is_operator = .True.
            case('-')
               is_operator = .True.
            case('*')
               is_operator = .True.
            case('/')
               is_operator = .True.
            case('^')
               is_operator = .True.
            case default
               is_operator = .False.
         end select

      end function is_operator


      !>
      !! @brief Function to decide whether a character is "*" or "/"
      !!
      !! ### Example
      !!~~~~~~~~~~~~~~.f90
      !! b = is_dot_operator("+")
      !! c = is_dot_operator("*")
      !!~~~~~~~~~~~~~~
      !! b will be .False. and c will be .True.
      !!
      !! @see is_line_operator, is_power_operator, is_separator, is_operator
      !!
      !! @author  Moritz Reichert
      !! @date    20.12.20
      function is_dot_operator(str_in)
         implicit none
         character(1),intent(in)  :: str_in          !< Input character
         logical                  :: is_dot_operator !< True in case of a * or / input string and False otherwise

         select case(str_in)
            case('*')
               is_dot_operator = .True.
            case('/')
               is_dot_operator = .True.
            case default
               is_dot_operator = .False.
         end select

      end function is_dot_operator


      !>
      !! @brief Function to decide whether a character is a "+" or "-"
      !!
      !! ### Example
      !!~~~~~~~~~~~~~~.f90
      !! b = is_line_operator("+")
      !! c = is_line_operator("*")
      !!~~~~~~~~~~~~~~
      !! b will be .True. and c will be .False.
      !!
      !! @see is_dot_operator, is_power_operator, is_separator, is_operator
      !!
      !! @author  Moritz Reichert
      !! @date    20.12.20
      function is_line_operator(str_in)
         implicit none
         character(1),intent(in)  :: str_in             !< Input character
         logical                  :: is_line_operator   !< True in case of a + or - input string and False otherwise

         select case(str_in)
            case('+')
               is_line_operator = .True.
            case('-')
               is_line_operator = .True.
            case default
               is_line_operator = .False.
         end select

      end function is_line_operator


      !>
      !! @brief Function to decide whether a character is a "^"
      !!
      !! ### Example
      !!~~~~~~~~~~~~~~.f90
      !! b = is_power_operator("^")
      !! c = is_power_operator("-")
      !!~~~~~~~~~~~~~~
      !! b will be .True. and c will be .False.
      !!
      !! @see is_dot_operator, is_line_operator, is_separator, is_operator
      !!
      !! @author  Moritz Reichert
      !! @date    20.12.20
      function is_power_operator(str_in)
         implicit none
         character(1),intent(in)  :: str_in              !< Input character
         logical                  :: is_power_operator   !< True in case of a ^ input string and False otherwise

         select case(str_in)
            case('^')
               is_power_operator = .True.
            case default
               is_power_operator = .False.
         end select

      end function is_power_operator


      !>
      !! @brief Function to decide whether a character can be associated with a number.
      !!
      !! Allowed characters are ".", "d", and "e"
      !! (for scientific format of a number).
      !!
      !! ### Example
      !!~~~~~~~~~~~~~~.f90
      !! b = is_separator(",")
      !! c = is_separator("-")
      !!~~~~~~~~~~~~~~
      !! b will be .True. and c will be .False.
      !!
      !! @see is_dot_operator, is_line_operator, is_separator, is_operator, is_power_operator
      !!
      !! @author  Moritz Reichert
      !! @date    20.12.20
      function is_separator(str_in)
         implicit none
         character(1),intent(in)  :: str_in          !< Input character
         logical                  :: is_separator    !< True in case of a ",", "d", and "e" input string and False otherwise

         select case(str_in)
            case('.')
               is_separator = .True.
            case('d')
               is_separator = .True.
            case('e')
               is_separator = .True.
            case default
               is_separator = .False.
         end select

      end function is_separator


      !>
      !! @brief Function to perform a simple mathematical operation.
      !!
      !! Depending on the input character, the function performs
      !! "+". "-", "*", "/", or "^". In case of a different character, an
      !! error is raised.
      !!
      !! ### Example
      !!~~~~~~~~~~~~~~.f90
      !! b = operation(10.3, 2.7, "+")
      !!~~~~~~~~~~~~~~
      !! b will be 13.
      !!
      !! @author  Moritz Reichert
      !! @date    20.12.20
      function operation(number_1,number_2,str_in)
         implicit none
         character(1),intent(in)  :: str_in               !< Input character (containing "+","-","*","/","^".)
         real(r_kind)             :: operation            !< Result of a certain operation of number_1 and number_2, determined by the input string.
         real(r_kind)             :: number_1,number_2    !< Input numbers on which the operation is performed

         select case(str_in)
            case('+')
               operation = number_1+number_2
            case('-')
               operation = number_1-number_2
            case('*')
               operation = number_1*number_2
            case('/')
               operation = number_1/number_2
            case('^')
               operation = number_1**number_2
            case default
               call raise_exception('Could not parse: '//trim(adjustl(complete_string))//&
                                    "."//NEW_LINE("A")//'Unknown operator "'//&
                                    trim(adjustl(str_in))//'".',"operation",&
                                    360003)
         end select

      end function operation


      !>
      !! @brief Function to perform a function operations.
      !!
      !! Depending on the input it will calculate,
      !! e.g., the square root (sqrt) or the sinus (sin).
      !! To add a function to the functionality of the parser,
      !! add the translation in this function here.
      !!
      !! ### Example
      !!~~~~~~~~~~~~~~.f90
      !! b = eval_function(4, "sqrt")
      !!~~~~~~~~~~~~~~
      !! b will be 2.
      !!
      !! @note New parsing functions can be defined here.
      !!
      !! @author  Moritz Reichert
      !! @date    20.12.20
      function eval_function(number,str_in)
         implicit none
         character(max_fname_len),intent(in)   :: str_in       !< Name of the function
         real(r_kind),intent(in)               :: number       !< input to the function
         real(r_kind)                          :: eval_function!< function(number), where function is defined by the input string.

         select case(trim(adjustl(str_in)))
            case('abs')
               eval_function =abs(number)
            case('sqrt')
               eval_function =sqrt(number)
            case('log')
               eval_function =log10(number)
            case('ln')
               eval_function =log(number)
            case('exp')
               eval_function =exp(number)
            case('sin')
               eval_function =sin(number)
            case('asin')
               eval_function =asin(number)
            case('cos')
               eval_function =cos(number)
            case('acos')
               eval_function =acos(number)
            case('tan')
               eval_function =tan(number)
            case('atan')
               eval_function =atan(number)
            case default
               call raise_exception('Could not parse: '//trim(adjustl(complete_string))//&
                                    "."//NEW_LINE("A")//'Function "'//trim(adjustl(str_in))//&
                                    '" not known.',"eval_function",360004)
         end select
      end function eval_function


      !>
      !! @brief Subroutine to check the parsing string for consistency and correct it.
      !!
      !! This subroutine will convert all characters to lower case characters
      !! (e.g., 4E-1 -> 4e-1) and will add a leading "0" in case the string starts
      !! with a mathematical operator (i.e., + or -, e.g., "-3+5" -> "0-3+5").
      !! Furthermore, it checks for "+-" (and "++", "--" and all combinations)
      !! and will correct it to "-". Additionally "**" is corrected to "^",
      !! in this way the  power can be written by both (**, and ^).
      !!
      !! ### Example
      !!~~~~~~~~~~~~~~.f90
      !! parsing_string = "-2+-3**2"
      !! call make_pars_str_consistent()
      !!~~~~~~~~~~~~~~
      !! After this, parsing_string will be "0-2-3^2"
      !!
      !! @returns A corrected string that can be parsed by other functions
      !!
      !! @author  Moritz Reichert
      !! @date    20.12.20
      subroutine make_pars_str_consistent
         implicit none
         integer           :: i,j              !< Loop variables
         integer           :: max_length_tmp   !< Length of the string without spaces
         logical           :: still_busy       !< Flag to check if the string still have to change

         max_length_tmp = len(trim(adjustl(parsing_string)))
         ! Convert the string to lower case letters to avoid problems
         do i = 1, max_length_tmp
           j = iachar(parsing_string(i:i))
           if (j>= iachar("A") .and. j<=iachar("Z") ) then
                parsing_string(i:i) = achar(iachar(parsing_string(i:i))+32)
           else
                cycle
           end if
         end do

         ! its a problem if the string starts with a minus or a plus
         if ((parsing_string(1:1) .eq. '-') .or. (parsing_string(1:1) .eq. '+')) then
            parsing_string = '0'//trim(adjustl(parsing_string))
         end if

         ! Loop through string to check for "++", "--", or "+-"
         still_busy = .True.
         do while (still_busy)
            still_busy = .False.
            do i=2,max_length_tmp
               ! +- is equal to -
               if (((parsing_string(i:i) .eq. '+') .and. (parsing_string(i-1:i-1) .eq. '-')) &
                  & .or.((parsing_string(i-1:i-1) .eq. '+') .and. (parsing_string(i:i) .eq. '-'))) then
                  parsing_string = parsing_string(1:i-2)//' -'//parsing_string(i+1:)
                  still_busy = .True.
               ! -- is equal to +
               else if ((parsing_string(i:i) .eq. '-') .and. (parsing_string(i-1:i-1) .eq. '-')) then
                  parsing_string = parsing_string(1:i-2)//' +'//parsing_string(i+1:)
                  still_busy = .True.
               ! ++ is equal to +
               else if ((parsing_string(i:i) .eq. '+') .and. (parsing_string(i-1:i-1) .eq. '+')) then
                  parsing_string = parsing_string(1:i-2)//' +'//parsing_string(i+1:)
                  still_busy = .True.
               ! ** is equal to ^
               else if ((parsing_string(i:i) .eq. '*') .and. (parsing_string(i-1:i-1) .eq. '*')) then
                  parsing_string = parsing_string(1:i-2)//' ^'//parsing_string(i+1:)
                  still_busy = .True.
               end if
            end do
         end do

      end subroutine make_pars_str_consistent


      !>
      !! @brief Function to define constants.
      !!
      !! This function translates names of constants
      !! to numbers (real). For example, "pi" is
      !! translated to 3.1415. Define new constants
      !! (names and values) in this function.
      !!
      !! ### Example
      !!~~~~~~~~~~~~~~.f90
      !! b = common_variables("pi")
      !!~~~~~~~~~~~~~~
      !! After this, will be 3.141592653589793
      !!
      !! @note New variables can be defined here.
      !!
      !! @author  Moritz Reichert
      !! @date    20.12.20
      function common_variables(inp_str)
         implicit none
         character(len=*),intent(in)  :: inp_str           !< Input character, defining the constant
         real(r_kind)                 :: common_variables  !< A value of a given constant

         ! Return common variables to replace them in the equation
         select case(trim(adjustl(inp_str)))
            case('pi')
               common_variables = 3.141592653589793
            case('kb')
               common_variables = 1.380649d-23
            case('e')
               common_variables = 2.718281828459045
            case(variable_name)
               common_variables = variable_value
            case default ! If it is not known variable
               common_variables = -1
         end select

      end function common_variables


      !>
      !! @brief Store mathematical operators (+-*/) and numbers in arrays.
      !!
      !! The result of this subroutine are two arrays, one containing all
      !! operators (e.g., for "3+2*2+1-1": ['+','*','+','-']) and one
      !! array containing all numbers (as real type, [3,2,2,1,1]).
      !!
      !! ### Example
      !!~~~~~~~~~~~~~~.f90
      !! parsing_string = "0-2-3^2"
      !! call store_digits_operator()
      !!~~~~~~~~~~~~~~
      !! After this, the array \ref operators will contain ("-","-","^"), and
      !! the array \ref numbers will contain (0, 2, 3, 2).
      !!
      !! @author  Moritz Reichert
      !! @date    20.12.20
      subroutine store_digits_operator
         implicit none
         integer                    :: i            !< Loop variable
         integer                    :: istat        !< status variable
         integer                    :: o_c_tmp      !< temporary operator count
         character(max_fname_len)   :: digit_tmp    !< temporary storage for numbers
         character(max_fname_len)   :: internal_pars!< internal parsing string
         real(r_kind)               :: tmp_storage  !< helper variable to replace constants

         ! Remove things from the parsing string and bring it
         ! to a shape that this subroutine can deal with
         call make_pars_str_consistent
         internal_pars = parsing_string

         ! Count the amount of operators
         operator_count = 0
         do i = 2,max_fname_len
            if ((is_operator(internal_pars(i:i))) &
               & .and. (.not. is_separator(internal_pars(i-1:i-1))) &
               & .and. (.not. is_operator(internal_pars(i-1:i-1)))) then
               operator_count = operator_count +1
            endif
         end do
         number_count = operator_count +1 ! amount of numbers

         ! deallocate if they are allocated
         if (allocated(operators)) then
            deallocate(operators)
            deallocate(numbers)
         end if

         allocate(operators(operator_count),numbers(number_count),stat=istat)
         if (istat /= 0) then
            call raise_exception('Could not parse: '//trim(adjustl(complete_string))//&
                                 "."//NEW_LINE("A")//&
                                 'Allocation of "operators" failed.',&
                                 "store_digits_operator",360001)
         end if
         o_c_tmp = 0
         digit_tmp = trim(adjustl(internal_pars(1:1))) ! It starts with a digit by defintion
         ! Count the amount of operators
         do i = 2,max_fname_len

            ! Store operators and numbers
            if ((is_operator(internal_pars(i:i)))&
               & .and. (.not. is_separator(internal_pars(i-1:i-1))) &
               & .and. (.not. is_operator(internal_pars(i-1:i-1)))) then
               o_c_tmp = o_c_tmp+1
               operators(o_c_tmp) = internal_pars(i:i)

               tmp_storage = common_variables(digit_tmp)
               if (tmp_storage .ne. -1) then ! It was a known variable name
                  numbers(o_c_tmp) = tmp_storage
               else
                  read( digit_tmp, *,iostat=istat) numbers(o_c_tmp)
                  if (istat /= 0) then ! Complain if something went wrong in conversion
                     call raise_exception('Could not parse: '//trim(adjustl(complete_string))//&
                                          "."//NEW_LINE("A")//&
                                          'Could not convert "'//trim(adjustl(digit_tmp))//'" to float.',&
                                          "store_digits_operator",360005)
                     end if
               end if
               digit_tmp = ''
            elseif (i .eq. max_fname_len-1) then ! Also store the last number
               o_c_tmp = o_c_tmp+1
               tmp_storage = common_variables(digit_tmp)
               if (tmp_storage .ne. -1) then ! It was a known variable name
                  numbers(o_c_tmp) = tmp_storage
               else ! otherwise read in the hopefully valid number
                  read( digit_tmp, *,iostat=istat) numbers(o_c_tmp)
                  if (istat /= 0) then ! Complain if something went wrong in conversion
                     call raise_exception('Could not parse: '//trim(adjustl(complete_string))//&
                                          "."//NEW_LINE("A")//&
                                          'Could not convert "'//trim(adjustl(digit_tmp))//'" to float.',&
                                          "store_digits_operator",360005)
                  end if
               end if
               digit_tmp = ''
            elseif ((.not. is_operator(internal_pars(i:i)))&
                   & .or. (is_separator(internal_pars(i-1:i-1))) &
                   & .or. (is_operator(internal_pars(i-1:i-1)))) then
               ! Append all characters to a digit (operators reset this string)
               digit_tmp = trim(adjustl(digit_tmp))//trim(adjustl(internal_pars(i:i)))
            elseif (internal_pars(i:i) .eq. ' ') then ! skip blanks in string
               cycle
            else
               call raise_exception('Could not parse: '//trim(adjustl(complete_string))//&
                                    "."//NEW_LINE("A")//&
                                    'Unknown character. ',&
                                    "store_digits_operator",360006)
            end if

         end do
      end subroutine store_digits_operator


      !>
      !! @brief Evaluate a simple string expression, without brackets and functions.
      !!
      !! This subroutine is able to evaluate simply expressions without
      !! brackets or functions. Only +- * / and ^ is allowed here. In addition
      !! a variable name or constant names may be given.
      !! For example "2+3-1" will store the result (4) in the first entry
      !! of the "numbers" array.
      !!
      !! ### Example
      !!~~~~~~~~~~~~~~.f90
      !! parsing_string = "0-2-3^2"
      !! call get_simple_result()
      !!~~~~~~~~~~~~~~
      !! After this, the array \ref numbers will have -11 on its first entry.
      !!
      !! @author  Moritz Reichert
      !! @date    20.12.20
      subroutine get_simple_result
         implicit none
         integer     :: i,j                            !< Loop variables
         integer     :: operator_count_tmp             !< Count variable
         logical     :: dot_operation                  !< still operations to do?
         logical     :: line_operation,power_operation !< still operations to do?
         real(r_kind):: helper_result                  !< Helper variable for evaluation

         if (.not. allocated(operators)) then
            call raise_exception('Could not parse: '//trim(adjustl(complete_string))//&
                                 "."//NEW_LINE("A")//&
                                 'Array "operators" was not allocated.',&
                                 "get_simple_result",360007)
         end if
         ! Calculate first exponents, then dot operations ( * / ) and then line operations (+ -)
         ! In German: "Punkt vor Strich ;)"

         operator_count_tmp = operator_count

         ! Power (^)
         power_operation = .True.
         do while(power_operation)
            power_operation = .False.
            innerloop_power: do i=1, operator_count_tmp
               ! Skip if all operators have been processed
               ! (the loop goes over the initial value of the variable)
               if (i .gt. operator_count_tmp) exit innerloop_power

               ! Found "^"
               if (is_power_operator(operators(i))) then
                  ! Evaluate the operation
                  helper_result = operation(numbers(i),numbers(i+1),operators(i))
                  numbers(i) = helper_result

                  ! shift all operators and numbers forward
                  do j=i,operator_count_tmp-1
                     operators(j) = operators(j+1)
                     numbers(j+1) = numbers(j+2)
                  end do
                  ! Still things to do
                  power_operation=.True.
                  ! Still things to do
                  operator_count_tmp = operator_count_tmp -1
               end if

            end do innerloop_power
         end do

         ! Dot operations (* /)
         dot_operation = .True.
         do while(dot_operation)
            dot_operation = .False.
            innerloop_dot: do i=1, operator_count_tmp
               ! Skip if all operators have been processed
               ! (the loop goes over the initial value of the variable)
               if (i .gt. operator_count_tmp) exit innerloop_dot
               ! Found * or /
               if (is_dot_operator(operators(i))) then
                  ! Evaluate the operation
                  helper_result = operation(numbers(i),numbers(i+1),operators(i))
                  numbers(i) = helper_result

                  ! shift all operators and numbers forward
                  do j=i,operator_count_tmp-1
                     operators(j) = operators(j+1)
                     numbers(j+1) = numbers(j+2)
                  end do
                  ! Still things to do
                  dot_operation=.True.
                  ! One operation has been performed
                  operator_count_tmp = operator_count_tmp -1
               end if

            end do innerloop_dot
         end do

         ! Line operations (+ -)
         line_operation = .True.
         do while(line_operation)
            line_operation = .False.
            innerloop_line: do i=1, operator_count_tmp
               ! Found + or -
               if (is_line_operator(operators(i))) then
                  ! Evaluate the operation
                  helper_result = operation(numbers(i),numbers(i+1),operators(i))
                  numbers(i) = helper_result

                  ! shift all operators and numbers forward
                  do j=i,operator_count_tmp-1
                     operators(j) = operators(j+1)
                     numbers(j+1) = numbers(j+2)
                  end do
                  ! Still things to do
                  line_operation=.True.
                  ! One operation has been performed
                  operator_count_tmp = operator_count_tmp -1
                  exit innerloop_line
               end if
            end do innerloop_line
         end do

      end subroutine get_simple_result


      !>
      !! @brief Evaluate a complex string expression.
      !!
      !! Evaluate expressions with brackets ( i.e., "(" and ")" ) and
      !! function names (i.e., "sin(5)" ).
      !!
      !! ### Example
      !!~~~~~~~~~~~~~~.f90
      !! parsing_string = "2+sin(pi)-5*(1+1)"
      !! call evaluate()
      !!~~~~~~~~~~~~~~
      !! After this, the array \ref numbers will have -8 on its first entry.
      !!
      !! @author  Moritz Reichert
      !! @date    20.12.20
      subroutine evaluate
         implicit none
         integer                           :: i !< Loop variable
         character(max_fname_len)          :: mod_expression
         character(1)                      :: cc   ! Current character in the loop
         integer                           :: start_block
         integer                           :: end_block,len_function_name
         character(max_fname_len)          :: string_res
         logical                           :: is_busy
         character(max_fname_len)          :: function_name,function_store

         mod_expression = complete_string
         is_busy=.True.

         ! Loop until nothing to do anymore
         do while (is_busy)
            is_busy = .False.
            function_name = ''
            inner_loop:do i=1,max_fname_len
               ! Store the current character to make it shorter
               cc = mod_expression(i:i)

               ! Ignore whitespaces
               if (cc .eq. ' ') then
                  cycle
               end if

               ! Keep track if a function is written in the expression
               if ((is_operator(cc) .or. is_digit(cc))) then
                  function_name = ''
               else
                  if (cc .ne. "(") then
                     function_name =trim(adjustl(function_name))//cc
                  end if
               end if

               ! Store functions and position of opening brackets
               if (mod_expression(i:i) .eq. '(') then
                  start_block = i+1
                  function_store=function_name
               end if

               ! The first closing bracket marks the innermost brackets
               if (mod_expression(i:i) .eq. ')') then
                  end_block = i-1
                  parsing_string = mod_expression(start_block:end_block)

                  ! Evaluate this sub-block that only has simple operators
                  call store_digits_operator
                  call get_simple_result

                  ! Should also a function get executed?
                  if (trim(adjustl(function_store)) .ne. '') then
                     numbers(1)=eval_function(numbers(1),function_store)
                     len_function_name = len(trim(adjustl(function_store))) ! get the length of the function name
                     start_block = start_block-len_function_name ! Also remove the function in the string
                  end if

                  write(string_res,*) numbers(1) ! The result is always stored in numbers(1) afterwards
                  mod_expression = mod_expression(1:start_block-2)//trim(adjustl(string_res))//mod_expression(end_block+2:max_fname_len)

                  ! write(*,*) trim(adjustl(mod_expression))
                  is_busy = .True.
                  exit inner_loop
               end if
            end do inner_loop
         end do

         ! Evaluate it a last time
         parsing_string = mod_expression
         call store_digits_operator
         call get_simple_result

         ! Now the final result is stored in numbers(1)
      end subroutine evaluate


      !>
      !! @brief Takes a string and evaluates the expression
      !!
      !! This function is accessible to the outside of the module.
      !! It serves as interface between other modules and this one.
      !! Also replaces \ref variable_name with \ref variable_value.
      !!
      !! ### Example
      !!~~~~~~~~~~~~~~.f90
      !! b = parse_string("-8+3*(2-1)-x",10)
      !!~~~~~~~~~~~~~~
      !! After this, b will be -15.
      !!
      !!
      !! @author  Moritz Reichert
      !! @date    20.12.20
      function parse_string(input_string,var_value)
         implicit none
         character(300),intent(in)  :: input_string !< String to parse
         real(r_kind),intent(in)    :: var_value    !< Value of the variable that will be stored in \ref variable_value
         real(r_kind)               :: parse_string !< Value of the evaluated expression

         INFO_ENTRY("PARSE_STRING")

         ! Store input in module variables
         variable_value = var_value
         complete_string = input_string

         ! Evaluate the expression
         call evaluate

         ! Store the result
         parse_string = numbers(1)

         ! deallocate the arrays if they are allocated
         if (allocated(operators)) then
            deallocate(operators)
            deallocate(numbers)
         end if

         INFO_EXIT("PARSE_STRING")
      end function parse_string


      !>
      !! @brief Finds a value of the variable for a given y-value
      !!
      !! This function solves an equation, i.e., finds the value of
      !! variable "x" for expressions as, for example, "3*x+2+x^2 = 10"
      !! returns 1.701562. This works with a newton raphson to find
      !! the root of f(x) - 10.
      !!
      !! ### Example
      !!~~~~~~~~~~~~~~.f90
      !! input_string  = "3*x+2+x^2"
      !! eq_value      = 10
      !! initial_guess = 0.0
      !! converged     = .False.
      !! result        = 0.0
      !! call find_start_value(input_string,eq_value,initial_guess,converged,result)
      !!~~~~~~~~~~~~~~
      !! After this, result will be 1.701562 and converged will be .True.
      !!
      !! @author  Moritz Reichert
      !! @date    20.12.20
      subroutine find_start_value(input_string,eq_value,initial_guess,converged,result)
         implicit none
         character(300),intent(in)  :: input_string       !< String to parse
         real(r_kind),intent(in)    :: eq_value           !< Value of the equation
         real(r_kind),intent(in)    :: initial_guess      !< Initial guess for the variable
         logical,intent(out)        :: converged          !< Flag that indicates a success
         real(r_kind),intent(out)   :: result             !< Time at which the expression is equal to eq_value.
         real(r_kind)               :: newguess           !< guess if the initial one failed
         integer                    :: i                  !< Loop variable
         integer,parameter          :: maxiter = 100      !< Maximum number of iterations
         real(r_kind)               :: y1,y2              !< values of f(x)-eq_value
         real(r_kind)               :: x1,x2,xnew         !< x values of the variable
         real(r_kind)               :: m,b                !< slope and intercept
         real(r_kind),parameter     :: diff = 1d-5        !< difference for the derivative
         real(r_kind),parameter     :: tol  = 1d-5        !< tolerance for the convergence
         integer                    :: try_count          !< count how often the initial guess was changed
         integer,parameter          :: max_try_count = 10 !< maximum changes of the initial guess

         INFO_ENTRY("FIND_START_VALUE")

         ! Set the first x-values
         newguess=initial_guess
         x1=newguess
         x2=x1+diff

         ! The outer loop is to ensure several initial guesses if NR fails
         converged = .False.
         outer_loop: do try_count=1, max_try_count

            ! Newton-Raphson
            do i=1, maxiter

               ! Calculate slope
               y1 = parse_string(input_string,x1)-eq_value
               y2 = parse_string(input_string,x2)-eq_value
               if (x1 .ne. x2) then
                  m  = (y1-y2)/(x1-x2)
               else
                  converged = .False.
                  exit outer_loop
               endif
               ! calculate new x value
               b  = y1-m*x1

               ! Exit if the slope is zero
               if (m .eq. 0) then
                  converged = .False.
                  ! Give a very large negative number
                  exit
               end if

               xnew = -b/m

               ! Exit if converged
               if (abs(x1-xnew) .lt. tol) then
                  converged = .True.
                  exit outer_loop
               end if

               ! Reshuffle variables if it is not converged
               x1 = xnew
               x2 = x1+diff
            end do

            ! Vary the initial guess
            newguess  = newguess+1d2
            x1        = newguess
            x2        = x1+diff
         end do outer_loop

      ! Store the result
      if (converged) result = xnew

      INFO_EXIT("FIND_START_VALUE")
      end subroutine find_start_value

end module parser_module
