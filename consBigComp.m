%constrains for big sized configuration while compression
function [c,ceq, gradc, gradceq] = consBigComp(x,delta,min)

ceq = [(x(42) - x(54))^2 + (x(43) - x(55))^2 + (x(44) - x(56))^2 - 0.3025, (x(18) - x(24))^2 + (x(19) - x(25))^2 + (x(20) - x(26))^2 - 0.3025, (-x(11) + x(8))^2 + (x(6) - x(9))^2 + (x(7) - x(10))^2 - 0.3025, (x(1) - x(12))^2 + (-x(13) + x(2))^2 + (-x(14) + x(3))^2 - 0.3025, (x(21) - x(63))^2 + (x(22) - x(64))^2 + (x(23) - x(65))^2 - 0.3025, (x(36) - x(51))^2 + (x(38) - x(53))^2 + (x(37) - x(52) + x(66))^2 - 0.3025, (x(34) - x(49))^2 + (x(35) - x(50))^2 + (delta + x(33) - x(48))^2 - 0.3025, (-x(57) + x(60))^2 + (-x(58) + x(61))^2 + (-x(59) + x(62) - x(67))^2 - 0.3025, (x(39) - x(45))^2 + (x(41) - x(47))^2 + (x(40) - x(46) - x(66))^2 - 0.3025, (-x(27) + x(30))^2 + (-x(28) + x(31))^2 + (-x(29) + x(32) + x(67))^2 - 0.3025, (-x(16) + 0.578524)^2 + (-x(17) + x(5))^2 + (-delta - x(15) + x(4))^2 - 0.3025, x(43)*(x(38) + x(41) + x(67)) - x(44)*(x(37) + x(40) - x(66)) + x(55)*(x(20) + x(29)) - x(56)*(x(19) + x(28) - x(66)), -x(42)*(x(38) + x(41) + x(67)) + x(44)*(x(36) + x(39)) - x(54)*(x(20) + x(29)) + x(56)*(x(18) + x(27)), x(42)*(x(37) + x(40) - x(66)) - x(43)*(x(36) + x(39)) + x(54)*(x(19) + x(28) - x(66)) - x(55)*(x(18) + x(27)), x(19)*(x(56) + 0.98557) - x(20)*(x(55) + x(66) + 0.578524) + x(25)*(x(11) + x(53)) - x(26)*(x(52) + x(10)), -x(18)*(x(56) + 0.98557) + x(20)*(x(54) + 0.389262) - x(24)*(x(11) + x(53)) + x(26)*(x(51) + x(9)), x(18)*(x(55) + x(66) + 0.578524) - x(19)*(x(54) + 0.389262) + x(24)*(x(52) + x(10)) - x(25)*(x(51) + x(9)), -x(11)*(x(25) + x(64)) + x(7)*(x(5) + 0.98557) - 2*0.578524*x(8) + x(10)*(x(26) + x(65)), x(11)*(delta + x(24) + x(63)) - x(6)*(x(5) + 0.98557) + x(8)*(x(4) + 0.939262) - x(9)*(x(26) + x(65)), 2*0.578524*x(6) - x(7)*(x(4) + 0.939262) + x(9)*(x(25) + x(64)) - x(10)*(delta + x(24) + x(63)), x(13)*(x(32) + x(62) - x(67)) - x(14)*(x(31) + x(61)) + x(2)*(x(5) - x(67) + 0.98557) - 2*0.578524*x(3), -x(1)*(x(5) - x(67) + 0.98557) - x(12)*(x(32) + x(62) - x(67)) + x(14)*(delta + x(30) + x(60)) + x(3)*(x(4) + 0.939262), 2*0.578524*x(1) + x(12)*(x(31) + x(61)) - x(13)*(delta + x(30) + x(60)) - x(2)*(x(4) + 0.939262), x(22)*(x(17) + x(50)) - x(23)*(x(16) + x(49) + x(66)) + x(64)*(x(11) + x(53)) - x(65)*(x(52) + x(10)), -x(21)*(x(17) + x(50)) + x(23)*(-delta + x(15) + x(48)) - x(63)*(x(11) + x(53)) + x(65)*(-delta + x(51) + x(9)), x(21)*(x(16) + x(49) + x(66)) - x(22)*(-delta + x(15) + x(48)) + x(63)*(x(52) + x(10)) - x(64)*(-delta + x(51) + x(9)), -x(38)*(x(34) + x(43) + 2*x(66)) + x(52)*(x(26) + x(65)) - x(53)*(x(25) + x(64)) + (x(35) + x(44))*(x(37) + x(66)), -x(36)*(x(35) + x(44)) + x(38)*(x(33) + x(42)) - x(51)*(x(26) + x(65)) + x(53)*(x(24) + x(63)), x(36)*(x(34) + x(43) + 2*x(66)) + x(51)*(x(25) + x(64)) - x(52)*(x(24) + x(63)) - (x(33) + x(42))*(x(37) + x(66)), x(34)*(x(38) + x(41)) - x(35)*(x(37) + x(40) - x(66)) + x(49)*(x(23) + x(59)) - x(50)*(x(22) + x(58) - x(66)), x(35)*(2*delta + x(36) + x(39)) - x(48)*(x(23) + x(59)) + x(50)*(2*delta + x(21) + x(57)) - (delta + x(33))*(x(38) + x(41)), -x(34)*(2*delta + x(36) + x(39)) + x(48)*(x(22) + x(58) - x(66)) - x(49)*(2*delta + x(21) + x(57)) + (delta + x(33))*(x(37) + x(40) - x(66)), x(58)*(x(17) + x(50) + 2*x(67)) + x(61)*(x(14) + x(47) + 2*x(67)) - x(62)*(x(13) + x(46)) - (x(16) + x(49))*(x(59) + x(67)), -x(57)*(x(17) + x(50) + 2*x(67)) - x(60)*(x(14) + x(47) + 2*x(67)) + x(62)*(-delta + x(12) + x(45)) + (x(59) + x(67))*(-delta + x(15) + x(48)), x(57)*(x(16) + x(49)) - x(58)*(-delta + x(15) + x(48)) + x(60)*(x(13) + x(46)) - x(61)*(-delta + x(12) + x(45)), x(40)*(x(35) + x(44) - x(67)) - x(41)*(x(34) + x(43) + 2*x(66)) - x(47)*(x(31) + x(61) + 2*x(66)) + (x(46) + x(66))*(x(32) + x(62) - x(67)), -x(39)*(x(35) + x(44) - x(67)) + x(41)*(x(33) + x(42)) - x(45)*(x(32) + x(62) - x(67)) + x(47)*(x(30) + x(60)), x(39)*(x(34) + x(43) + 2*x(66)) - x(40)*(x(33) + x(42)) + x(45)*(x(31) + x(61) + 2*x(66)) - (x(30) + x(60))*(x(46) + x(66)), x(28)*(x(56) + 0.98557) - x(29)*(x(55) + 0.578524) + x(31)*(x(14) + x(47) + 2*x(67)) - (x(13) + x(46))*(x(32) + x(67)), -x(27)*(x(56) + 0.98557) + x(29)*(x(54) + 0.389262) - x(30)*(x(14) + x(47) + 2*x(67)) + (x(12) + x(45))*(x(32) + x(67)), x(27)*(x(55) + 0.578524) - x(28)*(x(54) + 0.389262) + x(30)*(x(13) + x(46)) - x(31)*(x(12) + x(45)), x(16)*(x(23) + x(59)) - x(17)*(x(22) + x(58)) + 0.578524*x(3) - x(5)*(x(2) + x(7)) + 0.578524*x(8), x(17)*(2*delta + x(21) + x(57)) - x(4)*(x(3) + x(8)) + x(5)*(x(1) + x(6)) + (-delta - x(15))*(x(23) + x(59)), -0.578524*x(1) - x(16)*(2*delta + x(21) + x(57)) + x(4)*(x(2) + x(7)) - 0.578524*x(6) + (delta + x(15))*(x(22) + x(58))];
%c=-1*[(-x(1) + 0.939262)^2 + (-x(2) + 0.578524)^2 + (-x(3) - x(67) + 0.98557)^2 - min075, (-x(6) + 0.939262)^2 + (-x(7) + 0.578524)^2 + (-x(8) + 0.98557)^2 - min075, (-x(27) + 0.389262)^2 + (-x(28) + 0.578524)^2 + (-x(29) + 0.98557)^2 - min075, (-x(18) + 0.389262)^2 + (-x(19) + 0.578524)^2 + (-x(20) + 0.98557)^2 - min075, (x(36) - x(42))^2 + (x(37) - x(43))^2 + (x(38) - x(44))^2 - min075, (x(27) - x(54))^2 + (x(28) - x(55))^2 + (x(29) - x(56))^2 - min075, (x(18) - x(54))^2 + (x(20) - x(56))^2 + (x(19) - x(55) - x(66))^2 - min075, (-x(39) + x(42))^2 + (-x(40) + x(43) + x(66))^2 + (-x(41) + x(44) - x(67))^2 - min075, (x(11) - x(26))^2 + (-x(24) + x(9))^2 + (-x(25) + x(10))^2 - min075, (x(24) - x(51))^2 + (x(25) - x(52))^2 + (x(26) - x(53))^2 - min075, (x(11) - x(65))^2 + (-x(64) + x(10))^2 + (-delta - x(63) + x(9))^2 - min075, (x(4) - x(6))^2 + (x(5) - x(8))^2 + (-x(7) + 0.578524)^2 - min075, (x(1) - x(4))^2 + (x(2) - 0.578524)^2 + (x(3) - x(5))^2 - min075, (x(12) - x(30))^2 + (x(13) - x(31))^2 + (x(14) - x(32))^2 - min075, (x(15) - x(21))^2 + (x(16) - x(22))^2 + (x(17) - x(23))^2 - min075, (x(51) - x(63))^2 + (x(52) - x(64))^2 + (x(53) - x(65))^2 - min075, (-x(23) + x(50))^2 + (-delta - x(21) + x(48))^2 + (-x(22) + x(49) + x(66))^2 - min075, (x(33) - x(36))^2 + (x(34) - x(37))^2 + (x(35) - x(38))^2 - min075, (-x(49) + x(58))^2 + (-x(50) + x(59))^2 + (delta - x(48) + x(57))^2 - min075, (-x(33) + x(39))^2 + (-x(35) + x(41))^2 + (-x(34) + x(40) - x(66))^2 - min075, (x(13) - x(61))^2 + (-delta + x(12) - x(60))^2 + (x(14) - x(62) + x(67))^2 - min075, (x(45) - x(60))^2 + (x(46) - x(61))^2 + (x(47) - x(62) + x(67))^2 - min075, (x(15) - x(57))^2 + (x(16) - x(58))^2 + (x(17) - x(59))^2 - min075, (x(30) - x(45))^2 + (x(31) - x(46))^2 + (x(32) - x(47))^2 - min075];
c=-1*[(-x(1) + 0.939262)^2 + (-x(2) + 0.578524)^2 + (-x(3) - x(67) + 0.98557)^2-min , (-x(6) + 0.939262)^2 + (-x(7) + 0.578524)^2 + (-x(8) + 0.98557)^2-min , (-x(27) + 0.389262)^2 + (-x(28) + 0.578524)^2 + (-x(29) + 0.98557)^2-min , (-x(18) + 0.389262)^2 + (-x(19) + 0.578524)^2 + (-x(20) + 0.98557)^2 - min , (x(36) - x(42))^2 + (x(37) - x(43))^2 + (x(38) - x(44))^2 - min , (x(27) - x(54))^2 + (x(28) - x(55))^2 + (x(29) - x(56))^2 - min , (x(18) - x(54))^2 + (x(20) - x(56))^2 + (x(19) - x(55) - x(66))^2- min , (-x(39) + x(42))^2 + (-x(40) + x(43) + x(66))^2 + (-x(41) + x(44) - x(67))^2 - min , (x(11) - x(26))^2 + (-x(24) + x(9))^2 + (-x(25) + x(10))^2 - min , (x(24) - x(51))^2 + (x(25) - x(52))^2 + (x(26) - x(53))^2 - min , (x(11) - x(65))^2 + (-x(64) + x(10))^2 + (-delta - x(63) + x(9))^2 - min , (x(4) - x(6))^2 + (x(5) - x(8))^2 + (-x(7) + 0.578524)^2 - min , (x(1) - x(4))^2 + (x(2) - 0.578524)^2 + (x(3) - x(5))^2 - min , (x(12) - x(30))^2 + (x(13) - x(31))^2 + (x(14) - x(32))^2 - min , (x(15) - x(21))^2 + (x(16) - x(22))^2 + (x(17) - x(23))^2 - min , (x(51) - x(63))^2 + (x(52) - x(64))^2 + (x(53) - x(65))^2 - min , (-x(23) + x(50))^2 + (-delta - x(21) + x(48))^2 + (-x(22) + x(49) + x(66))^2 - min , (x(33) - x(36))^2 + (x(34) - x(37))^2 + (x(35) - x(38))^2 - min , (-x(49) + x(58))^2 + (-x(50) + x(59))^2 + (delta - x(48) + x(57))^2 - min , (-x(33) + x(39))^2 + (-x(35) + x(41))^2 + (-x(34) + x(40) - x(66))^2 - min , (x(13) - x(61))^2 + (-delta + x(12) - x(60))^2 + (x(14) - x(62) + x(67))^2 - min , (x(45) - x(60))^2 + (x(46) - x(61))^2 + (x(47) - x(62) + x(67))^2 - min , (x(15) - x(57))^2 + (x(16) - x(58))^2 + (x(17) - x(59))^2 - min , (x(30) - x(45))^2 + (x(31) - x(46))^2 + (x(32) - x(47))^2 - min ];
gradceq =[0, 0, 0, 2*x(1) - 2*x(12), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,        0, 0, 0, 0, -x(5) + x(67) - 0.98557, 1.15704800000000, 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x(5),         -0.578524000000000 ; 0, 0, 0, -2*x(13) + 2*x(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, x(5) - x(67) + 0.98557, 0, -x(4) - 0.939262, 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -x(5), 0, x(4) ; 0, 0, 0, -2*x(14) + 2*x(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, -1.15704800000000, x(4) + 0.939262, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.578524000000000, -x(4), 0 ; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2*delta - 2*x(15) + 2*x(4), 0, 0, 0,         0, 0, 0, 0, x(8), -x(7), 0, x(3), -x(2), 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -x(3) - x(8), x(2) + x(7) ; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2*x(17) + 2*x(5), 0, 0, 0, 0, 0, 0,         x(7), -x(6), 0, x(2), -x(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, -x(2) - x(7), x(1) + x(6), 0 ; 0, 0, 2*x(6) - 2*x(9), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         -x(5) - 0.98557, 1.15704800000000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x(5), -0.578524000000000 ; 0, 0, 2*x(7) - 2*x(10), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         x(5) + 0.98557, 0, -x(4) - 0.939262, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -x(5), 0, x(4) ; 0, 0, -2*x(11) + 2*x(8), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         -1.15704800000000, x(4) + 0.939262, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.578524000000000, -x(4), 0 ; 0, 0, -2*x(6) + 2*x(9), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x(26),         -x(25), 0, -x(26) - x(65), x(25) + x(64), 0, 0, 0, 0, x(65),         -x(64), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ; 0,0, -2*x(7) + 2*x(10), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -x(26), 0,         x(24), x(26) + x(65), 0, -delta - x(24) - x(63), 0, 0, 0, -x(65),         0, x(63), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ; 0, 0, 2*x(11) - 2*x(8), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x(25),         -x(24), 0, -x(25) - x(64), delta + x(24) + x(63), 0, 0, 0, 0,         x(64), -x(63), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, 0 ; 0, 0, 0, -2*x(1) + 2*x(12), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, -x(32) - x(62) + x(67), x(31) + x(61), 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, x(62), -x(61), 0, 0, 0, 0, x(32) + x(67), -x(31), 0,         0, 0 ; 0, 0, 0, 2*x(13) - 2*x(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, x(32) + x(62) - x(67), 0, -delta - x(30) - x(60), 0, 0, 0,         0, 0, 0, 0, 0, 0, -x(62), 0, x(60), 0, 0, 0, -x(32) - x(67), 0,         x(30), 0, 0, 0 ; 0, 0, 0, 2*x(14) - 2*x(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, -x(31) - x(61), delta + x(30) + x(60), 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, x(61), -x(60), 0, 0, 0, 0, x(31), -x(30), 0, 0, 0, 0 ; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2*delta + 2*x(15) - 2*x(4), 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x(23), -x(22), 0, 0, 0, 0, 0, 0, 0,         x(59) + x(67), -x(58), 0, 0, 0, 0, 0, 0, 0, -x(23) - x(59),         x(22) + x(58) ; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2*x(16) - 1.157048, 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, -x(23), 0, x(21), 0, 0, 0, 0, 0, 0,         -x(59) - x(67), 0, x(57), 0, 0, 0, 0, 0, 0, x(23) + x(59), 0,         -2*delta - x(21) - x(57) ; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2*x(17) - 2*x(5), 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, x(22), -x(21), 0, 0, 0, 0, 0, 0, 0, x(58),         -x(57), 0, 0, 0, 0, 0, 0, 0, -x(22) - x(58),         2*delta + x(21) + x(57), 0 ; 0, 2*x(18) - 2*x(24), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x(56), -x(55),         0, -x(56) - 0.98557, x(55) + x(66) + 0.578524, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ; 0, 2*x(19) - 2*x(25), 0, 0, 0, 0, 0, 0, 0, 0, 0, -x(56), 0, x(54),         x(56) + 0.98557, 0, -x(54) - 0.389262, 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ; 0, 2*x(20) - 2*x(26), 0, 0, 0, 0, 0, 0, 0, 0, 0, x(55), -x(54), 0,         -x(55) - x(66) - 0.578524, x(54) + 0.389262, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ; 0, 0, 0, 0, 2*x(21) - 2*x(63), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, -x(17) - x(50), x(16) + x(49) + x(66), 0, 0,         0, 0, x(50), -x(49), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x(17), -x(16) ; 0, 0, 0, 0, 2*x(22) - 2*x(64), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, x(17) + x(50), 0, delta - x(15) - x(48), 0, 0, 0,         -x(50), 0, x(48), 0, 0, 0, 0, 0, 0, 0, 0, 0, -x(17), 0,         delta + x(15) ; 0, 0, 0, 0, 2*x(23) - 2*x(65), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, -x(16) - x(49) - x(66), -delta + x(15) + x(48),         0, 0, 0, 0, x(49), -x(48), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x(16),         -delta - x(15), 0 ; 0, -2*x(18) + 2*x(24), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         -x(11) - x(53), x(52) + x(10), 0, x(11), -x(10), 0, 0, 0, 0, 0, 0, 0,         x(53), -x(52), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ; 0, -2*x(19) + 2*x(25), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         x(11) + x(53), 0, -x(51) - x(9), -x(11), 0, x(9), 0, 0, 0, 0, 0, 0,         -x(53), 0, x(51), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ; 0, -2*x(20) + 2*x(26), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         -x(52) - x(10), x(51) + x(9), 0, x(10), -x(9), 0, 0, 0, 0, 0, 0, 0,         x(52), -x(51), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ; 0, 0, 0, 0, 0, 0, 0, 0, 0, 2*x(27) - 2*x(30), 0, 0, x(56), -x(55),         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, -x(56) - 0.98557, x(55) + 0.578524, 0, 0, 0 ; 0, 0, 0, 0 ,0, 0, 0, 0, 0, 2*x(28) - 2*x(31), 0, -x(56), 0, x(54),         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, x(56) + 0.98557, 0, -x(54) - 0.389262, 0, 0, 0 ; 0, 0, 0, 0, 0, 0, 0, 0, 0, 2*x(29) - 2*x(32) - 2*x(67), 0, x(55),         -x(54), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, -x(55) - 0.578524, x(54) + 0.389262, 0, 0, 0, 0 ; 0, 0, 0, 0, 0, 0, 0, 0, 0, -2*x(27) + 2*x(30), 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, x(14), -x(13), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         x(47), -x(46) - x(66), 0, -x(14) - x(47) - 2*x(67), x(13) + x(46),         0, 0, 0 ; 0, 0, 0, 0, 0, 0, 0, 0, 0, -2*x(28) + 2*x(31), 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, -x(14), 0, x(12), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         -x(47), 0, x(45), x(14) + x(47) + 2*x(67), 0, -x(12) - x(45), 0, 0,         0 ; 0, 0, 0, 0, 0, 0, 0, 0, 0, -2*x(29) + 2*x(32) + 2*x(67), 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, x(13), -x(12), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, x(46) + x(66), -x(45), 0, -x(13) - x(46), x(12) + x(45),         0, 0, 0, 0 ; 0, 0, 0, 0, 0, 0, 2*delta + 2*x(33) - 2*x(48), 0, 0, 0, 0, 0, 0, 0,   0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x(38), -x(37) - x(66), 0,         -x(38) - x(41), x(37) + x(40) - x(66), 0, 0, 0, 0, x(41), -x(40),         0, 0, 0, 0, 0, 0 ; 0, 0, 0, 0, 0, 0, 2*x(34) - 2*x(49), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, -x(38), 0, x(36), x(38) + x(41), 0,         -2*delta - x(36) - x(39), 0, 0, 0, -x(41), 0, x(39), 0, 0, 0, 0, 0,         0 ; 0, 0, 0, 0, 0, 0, 2*x(35) - 2*x(50), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, x(37) + x(66), -x(36), 0,         -x(37) - x(40) + x(66), 2*delta + x(36) + x(39), 0, 0, 0, 0, x(40),         -x(39), 0, 0, 0, 0, 0, 0, 0 ; 0, 0, 0, 0, 0, 2*x(36) - 2*x(51), 0, 0, 0, 0, 0, 0, x(44), -x(43),         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -x(35) - x(44),         x(34) + x(43) + 2*x(66), 0, x(35), -x(34), 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, 0 ; 0, 0, 0, 0, 0, 2*x(37) - 2*x(52) + 2*x(66), 0, 0, 0, 0, 0, -x(44),         0, x(42), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x(35) + x(44), 0,         -x(33) - x(42), -x(35), 0, delta + x(33), 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, 0 ; 0, 0, 0, 0, 0, 2*x(38) - 2*x(53), 0, 0, 0, 0, 0, x(43), -x(42), 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -x(34) - x(43) - 2*x(66),         x(33) + x(42), 0, x(34), -delta - x(33), 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, 0 ; 0, 0, 0, 0, 0, 0, 0, 0, 2*x(39) - 2*x(45), 0, 0, 0, x(44), -x(43),         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x(35), -x(34), 0,         0, 0, 0, -x(35) - x(44) + x(67), x(34) + x(43) + 2*x(66), 0, 0, 0,         0, 0, 0 ; 0, 0, 0, 0, 0, 0, 0, 0, 2*x(40) - 2*x(46) - 2*x(66), 0, 0, -x(44),         0, x(42), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -x(35), 0,         delta + x(33), 0, 0, 0, x(35) + x(44) - x(67), 0, -x(33) - x(42),         0, 0, 0, 0, 0, 0 ; 0, 0, 0, 0, 0, 0, 0, 0, 2*x(41) - 2*x(47), 0, 0, x(43), -x(42), 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x(34), -delta - x(33),         0, 0, 0, 0, -x(34) - x(43) - 2*x(66), x(33) + x(42), 0, 0, 0, 0, 0,         0, 0 ; 2*x(42) - 2*x(54), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         -x(38) - x(41) - x(67), x(37) + x(40) - x(66), 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, x(38), -x(37) - x(66), 0, 0, 0, 0, 0, 0, 0,         x(41), -x(40), 0, 0, 0, 0, 0, 0 ; 2*x(43) - 2*x(55), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         x(38) + x(41) + x(67), 0, -x(36) - x(39), 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, -x(38), 0, x(36), 0, 0, 0, 0, 0, 0, -x(41), 0, x(39),         0, 0, 0, 0, 0, 0 ; 2*x(44) - 2*x(56), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         -x(37) - x(40) + x(66), x(36) + x(39), 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, x(37) + x(66), -x(36), 0, 0, 0, 0, 0, 0, 0, x(40),         -x(39), 0, 0, 0, 0, 0, 0, 0 ; 0, 0, 0, 0, 0, 0, 0, 0, -2*x(39) + 2*x(45), 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x(62), -x(61), 0,         -x(32) - x(62) + x(67), x(31) + x(61) + 2*x(66), 0, x(32) + x(67),         -x(31), 0, 0, 0 ; 0, 0, 0, 0, 0, 0, 0, 0, -2*x(40) + 2*x(46) + 2*x(66), 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -x(62), 0,         x(60), x(32) + x(62) - x(67), 0, -x(30) - x(60), -x(32) - x(67), 0,         x(30), 0, 0, 0 ; 0, 0, 0, 0, 0, 0, 0, 0, -2*x(41) + 2*x(47), 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x(61), -x(60), 0,         -x(31) - x(61) - 2*x(66), x(30) + x(60), 0, x(31), -x(30), 0, 0, 0,         0 ; 0, 0, 0, 0, 0, 0, -2*delta - 2*x(33) + 2*x(48), 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x(23), -x(22), 0, 0, 0, 0,         -x(23) - x(59), x(22) + x(58) - x(66), 0, x(59) + x(67), -x(58), 0,         0, 0, 0, 0, 0, 0, 0, 0 ; 0, 0, 0, 0, 0, 0, -2*x(34) + 2*x(49), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, -x(23), 0, x(21), 0, 0, 0, x(23) + x(59), 0,         -2*delta - x(21) - x(57), -x(59) - x(67), 0, x(57), 0, 0, 0, 0, 0,         0, 0, 0, 0 ; 0, 0, 0, 0, 0, 0, -2*x(35) + 2*x(50), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, x(22), -x(21), 0, 0, 0, 0,         -x(22) - x(58) + x(66), 2*delta + x(21) + x(57), 0, x(58), -x(57),         0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ; 0, 0, 0, 0, 0, -2*x(36) + 2*x(51), 0, 0, 0, 0, 0, 0, 0, 0, 0,         x(26), -x(25), 0, 0, 0, 0, 0, 0, 0, x(65), -x(64), 0,         -x(26) - x(65), x(25) + x(64), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0 ; 0, 0, 0, 0, 0, -2*x(37) + 2*x(52) - 2*x(66), 0, 0, 0, 0, 0, 0, 0,         0, -x(26), 0, x(24), 0, 0, 0, 0, 0, 0, -x(65), 0, x(63),         x(26) + x(65), 0, -x(24) - x(63), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, 0 ; 0, 0, 0, 0, 0, -2*x(38) + 2*x(53), 0, 0, 0, 0, 0, 0, 0, 0, x(25),         -x(24), 0, 0, 0, 0, 0, 0, 0, x(64), -x(63), 0, -x(25) - x(64),         x(24) + x(63), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ; -2*x(42) + 2*x(54), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         -x(20) - x(29), x(19) + x(28) - x(66), 0, x(20), -x(19), 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x(29),         -x(28), 0, 0, 0 ; -2*x(43) + 2*x(55), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x(20) + x(29), 0,         -x(18) - x(27), -x(20), 0, x(18), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -x(29), 0, x(27), 0, 0, 0 ; -2*x(44) + 2*x(56), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         -x(19) - x(28) + x(66), x(18) + x(27), 0, x(19), -x(18), 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x(28),         -x(27), 0, 0, 0, 0 ; 0, 0, 0, 0, 0, 0, 0, 2*x(57) - 2*x(60), 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x(50), -x(49), 0,         -x(17) - x(50) - 2*x(67), x(16) + x(49), 0, 0, 0, 0, 0, 0, 0,         x(17), -x(16) ; 0, 0, 0, 0, 0, 0, 0, 2*x(58) - 2*x(61), 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -x(50), 0, x(48),         x(17) + x(50) + 2*x(67), 0, delta - x(15) - x(48), 0, 0, 0, 0, 0,         0, -x(17), 0, delta + x(15) ; 0, 0, 0, 0, 0, 0, 0, 2*x(59) - 2*x(62) + 2*x(67), 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, x(49), -x(48), 0,         -x(16) - x(49), -delta + x(15) + x(48), 0, 0, 0, 0, 0, 0, 0, x(16),         -delta - x(15), 0 ; 0, 0, 0, 0, 0, 0, 0, -2*x(57) + 2*x(60), 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, x(14), -x(13), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         -x(14) - x(47) - 2*x(67), x(13) + x(46), 0, x(47), -x(46) - x(66),         0, 0, 0, 0, 0, 0 ; 0, 0, 0, 0, 0, 0, 0, -2*x(58) + 2*x(61), 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, 0, -x(14), 0, x(12), 0, 0, 0, 0, 0, 0, 0, 0, 0,         x(14) + x(47) + 2*x(67), 0, delta - x(12) - x(45), -x(47), 0,         x(45), 0, 0, 0, 0, 0, 0 ; 0, 0, 0, 0, 0, 0, 0, -2*x(59) + 2*x(62) - 2*x(67), 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, x(13), -x(12), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         -x(13) - x(46), -delta + x(12) + x(45), 0, x(46) + x(66), -x(45),         0, 0, 0, 0, 0, 0, 0 ; 0, 0, 0, 0, -2*x(21) + 2*x(63), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, x(11), -x(10), 0, 0, 0, 0, -x(11) - x(53), x(52) + x(10), 0,         x(53), -x(52), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ; 0, 0, 0, 0, -2*x(22) + 2*x(64), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         -x(11), 0, x(9), 0, 0, 0, x(11) + x(53), 0, delta - x(51) - x(9),         -x(53), 0, x(51), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ; 0, 0, 0, 0, -2*x(23) + 2*x(65), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         x(10), -x(9), 0, 0, 0, 0, -x(52) - x(10), -delta + x(51) + x(9), 0,         x(52), -x(51), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 ; 0, 0, 0, 0, 0, 2*x(37) - 2*x(52) + 2*x(66), 0, 0,         -2*x(40) + 2*x(46) + 2*x(66), 0, 0, x(44) + x(56), 0,         -x(42) - x(54), -x(20), 0, x(18), 0, 0, 0, 0, 0, 0, -x(23), 0,         x(21), x(35) - 2*x(38) + x(44), 0, -x(33) + 2*x(36) - x(42),         x(35) + x(50), 0, -delta - x(33) - x(48), 0, 0, 0,         x(32) - 2*x(41) - 2*x(47) + x(62) - x(67), 0,         -x(30) + 2*x(39) + 2*x(45) - x(60), 0, 0, 0, 0, 0, 0 ; 0, 0, 0, 0, 0, 0, 0, 2*x(59) - 2*x(62) + 2*x(67), 0,         -2*x(29) + 2*x(32) + 2*x(67), 0, x(43), -x(42), 0, 0, 0, 0, 0, 0,         0, -x(13) - x(2), x(1) + x(12), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         -x(16) - x(49) + 2*x(58) + 2*x(61),  -delta + x(15) + x(48) - 2*x(57) - 2*x(60), 0,  -x(40) - x(46) - x(66), x(39) + x(45), 0, -x(13) + 2*x(31) - x(46), x(12) - 2*x(30) + x(45), 0, 0, 0, 0  ];
gradc = [2*x(1) - 1.878524, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         2*x(1) - 2*x(4), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;2*x(2) - 1.157048, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         2*x(2) - 1.157048, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;2*x(3) + 2*x(67) - 1.97114, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         2*x(3) - 2*x(5), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2*x(4) - 2*x(6), -2*x(1) + 2*x(4),         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2*x(5) - 2*x(8), -2*x(3) + 2*x(5),         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;0, 2*x(6) - 1.878524, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2*x(4) + 2*x(6),         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;0, 2*x(7) - 1.157048, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2*x(7) - 1.157048,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;0, 2*x(8) - 1.97114, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2*x(5) + 2*x(8),         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;0, 0, 0, 0, 0, 0, 0, 0, -2*x(24) + 2*x(9), 0,         -2*delta - 2*x(63) + 2*x(9), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;0, 0, 0, 0, 0, 0, 0, 0, -2*x(25) + 2*x(10), 0, -2*x(64) + 2*x(10), 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;0, 0, 0, 0, 0, 0, 0, 0, 2*x(11) - 2*x(26), 0, 2*x(11) - 2*x(65), 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2*x(12) - 2*x(30), 0, 0, 0,         0, 0, 0, -2*delta + 2*x(12) - 2*x(60), 0, 0, 0;0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2*x(13) - 2*x(31), 0, 0, 0,         0, 0, 0, 2*x(13) - 2*x(61), 0, 0, 0;0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2*x(14) - 2*x(32), 0, 0, 0,         0, 0, 0, 2*x(14) - 2*x(62) + 2*x(67), 0, 0, 0;0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2*x(15) - 2*x(21), 0, 0,         0, 0, 0, 0, 0, 2*x(15) - 2*x(57), 0;0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2*x(16) - 2*x(22), 0, 0,         0, 0, 0, 0, 0, 2*x(16) - 2*x(58), 0;0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2*x(17) - 2*x(23), 0, 0,         0, 0, 0, 0, 0, 2*x(17) - 2*x(59), 0;0, 0, 0, 2*x(18) - 0.778524, 0, 0, 2*x(18) - 2*x(54), 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;0, 0, 0, 2*x(19) - 1.157048, 0, 0, 2*x(19) - 2*x(55) - 2*x(66), 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;0, 0, 0, 2*x(20) - 1.97114, 0, 0, 2*x(20) - 2*x(56), 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2*x(15) + 2*x(21), 0,         2*delta + 2*x(21) - 2*x(48), 0, 0, 0, 0, 0, 0, 0;0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2*x(16) + 2*x(22), 0,         2*x(22) - 2*x(49) - 2*x(66), 0, 0, 0, 0, 0, 0, 0;0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2*x(17) + 2*x(23), 0,         2*x(23) - 2*x(50), 0, 0, 0, 0, 0, 0, 0;0, 0, 0, 0, 0, 0, 0, 0, 2*x(24) - 2*x(9), 2*x(24) - 2*x(51), 0,0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;0, 0, 0, 0, 0, 0, 0, 0, 2*x(25) - 2*x(10), 2*x(25) - 2*x(52), 0, 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;0, 0, 0, 0, 0, 0, 0, 0, -2*x(11) + 2*x(26), 2*x(26) - 2*x(53), 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;0, 0, 2*x(27) - 0.778524, 0, 0, 2*x(27) - 2*x(54), 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;0, 0, 2*x(28) - 1.157048, 0, 0, 2*x(28) - 2*x(55), 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;0, 0, 2*x(29) - 1.97114, 0, 0, 2*x(29) - 2*x(56), 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2*x(12) + 2*x(30), 0, 0, 0,         0, 0, 0, 0, 0, 0, 2*x(30) - 2*x(45);0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2*x(13) + 2*x(31), 0, 0, 0,         0, 0, 0, 0, 0, 0, 2*x(31) - 2*x(46);0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2*x(14) + 2*x(32), 0, 0, 0,         0, 0, 0, 0, 0, 0, 2*x(32) - 2*x(47);0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         2*x(33) - 2*x(36), 0, 2*x(33) - 2*x(39), 0, 0, 0, 0;0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         2*x(34) - 2*x(37), 0, 2*x(34) - 2*x(40) + 2*x(66), 0, 0, 0, 0;0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         2*x(35) - 2*x(38), 0, 2*x(35) - 2*x(41), 0, 0, 0, 0;0, 0, 0, 0, 2*x(36) - 2*x(42), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         -2*x(33) + 2*x(36), 0, 0, 0, 0, 0, 0;0, 0, 0, 0, 2*x(37) - 2*x(43), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         -2*x(34) + 2*x(37), 0, 0, 0, 0, 0, 0;0, 0, 0, 0, 2*x(38) - 2*x(44), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         -2*x(35) + 2*x(38), 0, 0, 0, 0, 0, 0;0, 0, 0, 0, 0, 0, 0, 2*x(39) - 2*x(42), 0, 0, 0, 0, 0, 0, 0, 0, 0,         0, 0, -2*x(33) + 2*x(39), 0, 0, 0, 0;0, 0, 0, 0, 0, 0, 0, 2*x(40) - 2*x(43) - 2*x(66), 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, -2*x(34) + 2*x(40) - 2*x(66), 0, 0, 0, 0;0, 0, 0, 0, 0, 0, 0, 2*x(41) - 2*x(44) + 2*x(67), 0, 0, 0, 0, 0, 0,         0, 0, 0, 0, 0, -2*x(35) + 2*x(41), 0, 0, 0, 0;0, 0, 0, 0, -2*x(36) + 2*x(42), 0, 0, -2*x(39) + 2*x(42), 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;0, 0, 0, 0, -2*x(37) + 2*x(43), 0, 0, -2*x(40) + 2*x(43) + 2*x(66),         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;0, 0, 0, 0, -2*x(38) + 2*x(44), 0, 0, -2*x(41) + 2*x(44) - 2*x(67),         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         2*x(45) - 2*x(60), 0, -2*x(30) + 2*x(45);0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         2*x(46) - 2*x(61), 0, -2*x(31) + 2*x(46);0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         2*x(47) - 2*x(62) + 2*x(67), 0, -2*x(32) + 2*x(47);0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         -2*delta - 2*x(21) + 2*x(48), 0, -2*delta + 2*x(48) - 2*x(57), 0,         0, 0, 0, 0;0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         -2*x(22) + 2*x(49) + 2*x(66), 0, 2*x(49) - 2*x(58), 0, 0, 0, 0, 0;0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2*x(23) + 2*x(50),         0, 2*x(50) - 2*x(59), 0, 0, 0, 0, 0;0, 0, 0, 0, 0, 0, 0, 0, 0, -2*x(24) + 2*x(51), 0, 0, 0, 0, 0,         2*x(51) - 2*x(63), 0, 0, 0, 0, 0, 0, 0, 0;0, 0, 0, 0, 0, 0, 0, 0, 0, -2*x(25) + 2*x(52), 0, 0, 0, 0, 0,         2*x(52) - 2*x(64), 0, 0, 0, 0, 0, 0, 0, 0;0, 0, 0, 0, 0, 0, 0, 0, 0, -2*x(26) + 2*x(53), 0, 0, 0, 0, 0,         2*x(53) - 2*x(65), 0, 0, 0, 0, 0, 0, 0, 0;0, 0, 0, 0, 0, -2*x(27) + 2*x(54), -2*x(18) + 2*x(54), 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;0, 0, 0, 0, 0, -2*x(28) + 2*x(55), -2*x(19) + 2*x(55) + 2*x(66), 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;0, 0, 0, 0, 0, -2*x(29) + 2*x(56), -2*x(20) + 2*x(56), 0, 0, 0, 0,         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         2*delta - 2*x(48) + 2*x(57), 0, 0, 0, -2*x(15) + 2*x(57), 0;0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         -2*x(49) + 2*x(58), 0, 0, 0, -2*x(16) + 2*x(58), 0;0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         -2*x(50) + 2*x(59), 0, 0, 0, -2*x(17) + 2*x(59), 0;0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         2*delta - 2*x(12) + 2*x(60), -2*x(45) + 2*x(60), 0, 0;0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         -2*x(13) + 2*x(61), -2*x(46) + 2*x(61), 0, 0;0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         -2*x(14) + 2*x(62) - 2*x(67), -2*x(47) + 2*x(62) - 2*x(67), 0, 0;0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2*delta + 2*x(63) - 2*x(9), 0, 0, 0,         0, -2*x(51) + 2*x(63), 0, 0, 0, 0, 0, 0, 0, 0;0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2*x(64) - 2*x(10), 0, 0, 0, 0,         -2*x(52) + 2*x(64), 0, 0, 0, 0, 0, 0, 0, 0;0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2*x(11) + 2*x(65), 0, 0, 0, 0,         -2*x(53) + 2*x(65), 0, 0, 0, 0, 0, 0, 0, 0;0, 0, 0, 0, 0, 0, -2*x(19) + 2*x(55) + 2*x(66),         -2*x(40) + 2*x(43) + 2*x(66), 0, 0, 0, 0, 0, 0, 0, 0,         -2*x(22) + 2*x(49) + 2*x(66), 0, 0, 2*x(34) - 2*x(40) + 2*x(66), 0,         0, 0, 0;2*x(3) + 2*x(67) - 1.97114, 0, 0, 0, 0, 0, 0,         2*x(41) - 2*x(44) + 2*x(67), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,         2*x(14) - 2*x(62) + 2*x(67), 2*x(47) - 2*x(62) + 2*x(67), 0, 0 ];
end
