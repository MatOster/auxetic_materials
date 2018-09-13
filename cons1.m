%constraints of medium sized system dependend on periodicity in x-direciton
function [c,ceq] = cons1(x,delta)
%delta =1.0001;
ceq = [(x(42) - x(54))^2 + (x(43) - x(55))^2 + (x(44) - x(56))^2 - 0.3025, (x(18) - x(24))^2 + (x(19) - x(25))^2 + (x(20) - x(26))^2 - 0.3025, (-x(11) + x(8))^2 + (x(6) - x(9))^2 + (x(7) - x(10))^2 - 0.3025, (x(1) - x(12))^2 + (-x(13) + x(2))^2 + (-x(14) + x(3))^2 - 0.3025, (x(21) - x(63))^2 + (x(22) - x(64))^2 + (x(23) - x(65))^2 - 0.3025, (x(36) - x(51))^2 + (x(38) - x(53))^2 + (x(37) - x(52) + x(66))^2 - 0.3025, (x(34) - x(49))^2 + (x(35) - x(50))^2 + (delta + x(33) - x(48))^2 - 0.3025, (-x(57) + x(60))^2 + (-x(58) + x(61))^2 + (-x(59) + x(62) - x(67))^2 - 0.3025, (x(39) - x(45))^2 + (x(41) - x(47))^2 + (x(40) - x(46) - x(66))^2 - 0.3025, (-x(27) + x(30))^2 + (-x(28) + x(31))^2 + (-x(29) + x(32) + x(67))^2 - 0.3025, (-x(16) + 0.5)^2 + (-x(17) + x(5))^2 + (-delta - x(15) + x(4))^2 - 0.3025, x(43)*(x(38) + x(41) + x(67)) - x(44)*(x(37) + x(40) - x(66)) + x(55)*(x(20) + x(29)) - x(56)*(x(19) + x(28) - x(66)), -x(42)*(x(38) + x(41) + x(67)) + x(44)*(x(36) + x(39)) - x(54)*(x(20) + x(29)) + x(56)*(x(18) + x(27)), x(42)*(x(37) + x(40) - x(66)) - x(43)*(x(36) + x(39)) + x(54)*(x(19) + x(28) - x(66)) - x(55)*(x(18) + x(27)), x(19)*(x(56) + 0.75) - x(20)*(x(55) + x(66) + 0.5) + x(25)*(x(11) + x(53)) - x(26)*(x(52) + x(10)), -x(18)*(x(56) + 0.75) + x(20)*(x(54) + 0.35) - x(24)*(x(11) + x(53)) + x(26)*(x(51) + x(9)), x(18)*(x(55) + x(66) + 0.5) - x(19)*(x(54) + 0.35) + x(24)*(x(52) + x(10)) - x(25)*(x(51) + x(9)), -x(11)*(x(25) + x(64)) + x(7)*(x(5) + 0.75) - x(8) + x(10)*(x(26) + x(65)), x(11)*(delta + x(24) + x(63)) - x(6)*(x(5) + 0.75) + x(8)*(x(4) + 0.9) - x(9)*(x(26) + x(65)), x(6) - x(7)*(x(4) + 0.9) + x(9)*(x(25) + x(64)) - x(10)*(delta + x(24) + x(63)), x(13)*(x(32) + x(62) - x(67)) - x(14)*(x(31) + x(61)) + x(2)*(x(5) - x(67) + 0.75) - x(3), -x(1)*(x(5) - x(67) + 0.75) - x(12)*(x(32) + x(62) - x(67)) + x(14)*(delta + x(30) + x(60)) + x(3)*(x(4) + 0.9), x(1) + x(12)*(x(31) + x(61)) - x(13)*(delta + x(30) + x(60)) - x(2)*(x(4) + 0.9), x(22)*(x(17) + x(50)) - x(23)*(x(16) + x(49) + x(66)) + x(64)*(x(11) + x(53)) - x(65)*(x(52) + x(10)), -x(21)*(x(17) + x(50)) + x(23)*(-delta + x(15) + x(48)) - x(63)*(x(11) + x(53)) + x(65)*(-delta + x(51) + x(9)), x(21)*(x(16) + x(49) + x(66)) - x(22)*(-delta + x(15) + x(48)) + x(63)*(x(52) + x(10)) - x(64)*(-delta + x(51) + x(9)), -x(38)*(x(34) + x(43) + 2*x(66)) + x(52)*(x(26) + x(65)) - x(53)*(x(25) + x(64)) + (x(35) + x(44))*(x(37) + x(66)), -x(36)*(x(35) + x(44)) + x(38)*(x(33) + x(42)) - x(51)*(x(26) + x(65)) + x(53)*(x(24) + x(63)), x(36)*(x(34) + x(43) + 2*x(66)) + x(51)*(x(25) + x(64)) - x(52)*(x(24) + x(63)) - (x(33) + x(42))*(x(37) + x(66)), x(34)*(x(38) + x(41)) - x(35)*(x(37) + x(40) - x(66)) + x(49)*(x(23) + x(59)) - x(50)*(x(22) + x(58) - x(66)), x(35)*(2*delta + x(36) + x(39)) - x(48)*(x(23) + x(59)) + x(50)*(2*delta + x(21) + x(57)) - (delta + x(33))*(x(38) + x(41)), -x(34)*(2*delta + x(36) + x(39)) + x(48)*(x(22) + x(58) - x(66)) - x(49)*(2*delta + x(21) + x(57)) + (delta + x(33))*(x(37) + x(40) - x(66)), x(58)*(x(17) + x(50) + 2*x(67)) + x(61)*(x(14) + x(47) + 2*x(67)) - x(62)*(x(13) + x(46)) - (x(16) + x(49))*(x(59) + x(67)), -x(57)*(x(17) + x(50) + 2*x(67)) - x(60)*(x(14) + x(47) + 2*x(67)) + x(62)*(-delta + x(12) + x(45)) + (x(59) + x(67))*(-delta + x(15) + x(48)), x(57)*(x(16) + x(49)) - x(58)*(-delta + x(15) + x(48)) + x(60)*(x(13) + x(46)) - x(61)*(-delta + x(12) + x(45)), x(40)*(x(35) + x(44) - x(67)) - x(41)*(x(34) + x(43) + 2*x(66)) - x(47)*(x(31) + x(61) + 2*x(66)) + (x(46) + x(66))*(x(32) + x(62) - x(67)), -x(39)*(x(35) + x(44) - x(67)) + x(41)*(x(33) + x(42)) - x(45)*(x(32) + x(62) - x(67)) + x(47)*(x(30) + x(60)), x(39)*(x(34) + x(43) + 2*x(66)) - x(40)*(x(33) + x(42)) + x(45)*(x(31) + x(61) + 2*x(66)) - (x(30) + x(60))*(x(46) + x(66)), x(28)*(x(56) + 0.75) - x(29)*(x(55) + 0.5) + x(31)*(x(14) + x(47) + 2*x(67)) - (x(13) + x(46))*(x(32) + x(67)), -x(27)*(x(56) + 0.75) + x(29)*(x(54) + 0.35) - x(30)*(x(14) + x(47) + 2*x(67)) + (x(12) + x(45))*(x(32) + x(67)), x(27)*(x(55) + 0.5) - x(28)*(x(54) + 0.35) + x(30)*(x(13) + x(46)) - x(31)*(x(12) + x(45)), x(16)*(x(23) + x(59)) - x(17)*(x(22) + x(58)) + 0.5*x(3) - x(5)*(x(2) + x(7)) + 0.5*x(8), x(17)*(2*delta + x(21) + x(57)) - x(4)*(x(3) + x(8)) + x(5)*(x(1) + x(6)) + (-delta - x(15))*(x(23) + x(59)), -0.5*x(1) - x(16)*(2*delta + x(21) + x(57)) + x(4)*(x(2) + x(7)) - 0.5*x(6) + (delta + x(15))*(x(22) + x(58))];
c=-1*[(-x(1) + 0.9)^2 + (-x(2) + 0.5)^2 + (-x(3) - x(67) + 0.75)^2 - 0.1075, (-x(6) + 0.9)^2 + (-x(7) + 0.5)^2 + (-x(8) + 0.75)^2 - 0.1075, (-x(27) + 0.35)^2 + (-x(28) + 0.5)^2 + (-x(29) + 0.75)^2 - 0.1075, (-x(18) + 0.35)^2 + (-x(19) + 0.5)^2 + (-x(20) + 0.75)^2 - 0.1075, (x(36) - x(42))^2 + (x(37) - x(43))^2 + (x(38) - x(44))^2 - 0.1075, (x(27) - x(54))^2 + (x(28) - x(55))^2 + (x(29) - x(56))^2 - 0.1075, (x(18) - x(54))^2 + (x(20) - x(56))^2 + (x(19) - x(55) - x(66))^2 - 0.1075, (-x(39) + x(42))^2 + (-x(40) + x(43) + x(66))^2 + (-x(41) + x(44) - x(67))^2 - 0.1075, (x(11) - x(26))^2 + (-x(24) + x(9))^2 + (-x(25) + x(10))^2 - 0.1075, (x(24) - x(51))^2 + (x(25) - x(52))^2 + (x(26) - x(53))^2 - 0.1075, (x(11) - x(65))^2 + (-x(64) + x(10))^2 + (-delta - x(63) + x(9))^2 - 0.1075, (x(4) - x(6))^2 + (x(5) - x(8))^2 + (-x(7) + 0.5)^2 - 0.1075, (x(1) - x(4))^2 + (x(2) - 0.5)^2 + (x(3) - x(5))^2 - 0.1075, (x(12) - x(30))^2 + (x(13) - x(31))^2 + (x(14) - x(32))^2 - 0.1075, (x(15) - x(21))^2 + (x(16) - x(22))^2 + (x(17) - x(23))^2 - 0.1075, (x(51) - x(63))^2 + (x(52) - x(64))^2 + (x(53) - x(65))^2 - 0.1075, (-x(23) + x(50))^2 + (-delta - x(21) + x(48))^2 + (-x(22) + x(49) + x(66))^2 - 0.1075, (x(33) - x(36))^2 + (x(34) - x(37))^2 + (x(35) - x(38))^2 - 0.1075, (-x(49) + x(58))^2 + (-x(50) + x(59))^2 + (delta - x(48) + x(57))^2 - 0.1075, (-x(33) + x(39))^2 + (-x(35) + x(41))^2 + (-x(34) + x(40) - x(66))^2 - 0.1075, (x(13) - x(61))^2 + (-delta + x(12) - x(60))^2 + (x(14) - x(62) + x(67))^2 - 0.1075, (x(45) - x(60))^2 + (x(46) - x(61))^2 + (x(47) - x(62) + x(67))^2 - 0.1075, (x(15) - x(57))^2 + (x(16) - x(58))^2 + (x(17) - x(59))^2 - 0.1075, (x(30) - x(45))^2 + (x(31) - x(46))^2 + (x(32) - x(47))^2 - 0.1075];
%c=[];
end
