function out1 = optimizedDynamics_SixBar_19_08_28_15_53_23(in1,in2,in3,in4)
%OPTIMIZEDDYNAMICS_SIXBAR_19_08_28_15_53_23
%    OUT1 = OPTIMIZEDDYNAMICS_SIXBAR_19_08_28_15_53_23(IN1,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.2.
%    28-Aug-2019 16:08:42

L1 = in4(1,:);
L2 = in4(2,:);
L3 = in4(3,:);
L4 = in4(4,:);
L5 = in4(5,:);
L6 = in4(6,:);
RL1 = in3(1,:);
RL2 = in3(2,:);
RL3 = in3(3,:);
RL4 = in3(4,:);
RL5 = in3(5,:);
RL6 = in3(6,:);
RL7 = in3(7,:);
RL8 = in3(8,:);
RL9 = in3(9,:);
RL10 = in3(10,:);
RL11 = in3(11,:);
RL12 = in3(12,:);
RL13 = in3(13,:);
RL14 = in3(14,:);
RL15 = in3(15,:);
RL16 = in3(16,:);
RL17 = in3(17,:);
RL18 = in3(18,:);
RL19 = in3(19,:);
RL20 = in3(20,:);
RL21 = in3(21,:);
RL22 = in3(22,:);
RL23 = in3(23,:);
RL24 = in3(24,:);
p1 = in1(1,:);
p2 = in1(2,:);
p3 = in1(3,:);
p4 = in1(4,:);
p5 = in1(5,:);
p6 = in1(6,:);
p7 = in1(7,:);
p8 = in1(8,:);
p9 = in1(9,:);
p10 = in1(10,:);
p11 = in1(11,:);
p12 = in1(12,:);
p13 = in1(13,:);
p14 = in1(14,:);
p15 = in1(15,:);
p16 = in1(16,:);
p17 = in1(17,:);
p18 = in1(18,:);
p19 = in1(19,:);
p20 = in1(20,:);
p21 = in1(21,:);
p22 = in1(22,:);
p23 = in1(23,:);
p24 = in1(24,:);
p25 = in1(25,:);
p26 = in1(26,:);
p27 = in1(27,:);
p28 = in1(28,:);
p29 = in1(29,:);
p30 = in1(30,:);
p31 = in1(31,:);
p32 = in1(32,:);
p33 = in1(33,:);
p34 = in1(34,:);
p35 = in1(35,:);
p36 = in1(36,:);
pDOT1 = in2(1,:);
pDOT2 = in2(2,:);
pDOT3 = in2(3,:);
pDOT4 = in2(4,:);
pDOT5 = in2(5,:);
pDOT6 = in2(6,:);
pDOT7 = in2(7,:);
pDOT8 = in2(8,:);
pDOT9 = in2(9,:);
pDOT10 = in2(10,:);
pDOT11 = in2(11,:);
pDOT12 = in2(12,:);
pDOT13 = in2(13,:);
pDOT14 = in2(14,:);
pDOT15 = in2(15,:);
pDOT16 = in2(16,:);
pDOT17 = in2(17,:);
pDOT18 = in2(18,:);
pDOT19 = in2(19,:);
pDOT20 = in2(20,:);
pDOT21 = in2(21,:);
pDOT22 = in2(22,:);
pDOT23 = in2(23,:);
pDOT24 = in2(24,:);
pDOT25 = in2(25,:);
pDOT26 = in2(26,:);
pDOT27 = in2(27,:);
pDOT28 = in2(28,:);
pDOT29 = in2(29,:);
pDOT30 = in2(30,:);
pDOT31 = in2(31,:);
pDOT32 = in2(32,:);
pDOT33 = in2(33,:);
pDOT34 = in2(34,:);
pDOT35 = in2(35,:);
pDOT36 = in2(36,:);
t2 = p3.^2;
t3 = p1-p19;
t4 = abs(t3);
t7 = p2-p20;
t5 = abs(t7);
t8 = p3-p21;
t6 = abs(t8);
t9 = t3.^2;
t10 = t7.^2;
t11 = t8.^2;
t12 = t9+t10+t11;
t13 = sqrt(t12);
t14 = RL5-t13;
t15 = p1-p22;
t16 = abs(t15);
t19 = p2-p23;
t17 = abs(t19);
t20 = p3-p24;
t18 = abs(t20);
t21 = t15.^2;
t22 = t19.^2;
t23 = t20.^2;
t24 = t21+t22+t23;
t25 = sqrt(t24);
t26 = RL6-t25;
t27 = p1-p28;
t28 = abs(t27);
t31 = p2-p29;
t29 = abs(t31);
t32 = p3-p30;
t30 = abs(t32);
t33 = t27.^2;
t34 = t31.^2;
t35 = t32.^2;
t36 = t33+t34+t35;
t37 = sqrt(t36);
t38 = RL7-t37;
t39 = p1-p34;
t40 = abs(t39);
t43 = p2-p35;
t41 = abs(t43);
t44 = p3-p36;
t42 = abs(t44);
t45 = t39.^2;
t46 = t43.^2;
t47 = t44.^2;
t48 = t45+t46+t47;
t49 = sqrt(t48);
t50 = RL8-t49;
t51 = p1-p4;
t52 = p1.^2;
t53 = p2.^2;
t54 = p4.^2;
t55 = p5.^2;
t56 = p6.^2;
t127 = p1.*p4.*2.0;
t128 = p2.*p5.*2.0;
t129 = p3.*p6.*2.0;
t57 = t2+t52+t53+t54+t55+t56-t127-t128-t129;
t58 = 1.0./t57;
t59 = p1.*2.0;
t60 = p4.*2.0;
t61 = t59-t60;
t62 = p2.*2.0;
t63 = p5.*2.0;
t64 = t62-t63;
t65 = p3.*2.0;
t66 = p6.*2.0;
t67 = t65-t66;
t68 = pDOT1.*2.0;
t69 = pDOT4.*2.0;
t70 = t68-t69;
t71 = pDOT2.*2.0;
t72 = pDOT5.*2.0;
t73 = t71-t72;
t74 = pDOT3.*2.0;
t75 = pDOT6.*2.0;
t76 = t74-t75;
t77 = p2-p5;
t78 = p3-p6;
t79 = p4-p13;
t80 = abs(t79);
t83 = p5-p14;
t81 = abs(t83);
t84 = p6-p15;
t82 = abs(t84);
t85 = t79.^2;
t86 = t83.^2;
t87 = t84.^2;
t88 = t85+t86+t87;
t89 = sqrt(t88);
t90 = RL2-t89;
t91 = p4-p16;
t92 = abs(t91);
t95 = p5-p17;
t93 = abs(t95);
t96 = p6-p18;
t94 = abs(t96);
t97 = t91.^2;
t98 = t95.^2;
t99 = t96.^2;
t100 = t97+t98+t99;
t101 = sqrt(t100);
t102 = RL1-t101;
t103 = p4-p28;
t104 = abs(t103);
t107 = p5-p29;
t105 = abs(t107);
t108 = p6-p30;
t106 = abs(t108);
t109 = t103.^2;
t110 = t107.^2;
t111 = t108.^2;
t112 = t109+t110+t111;
t113 = sqrt(t112);
t114 = RL3-t113;
t115 = p4-p34;
t116 = abs(t115);
t119 = p5-p35;
t117 = abs(t119);
t120 = p6-p36;
t118 = abs(t120);
t121 = t115.^2;
t122 = t119.^2;
t123 = t120.^2;
t124 = t121+t122+t123;
t125 = sqrt(t124);
t126 = RL4-t125;
t130 = t56+1.0e-4;
t131 = 1.0./sqrt(t130);
t132 = (p6.*t131)./2.0;
t133 = t132-1.0./2.0;
t134 = t80.^2;
t135 = t134.*2.0;
t136 = t81.^2;
t137 = t136.*2.0;
t138 = t82.^2;
t139 = t138.*2.0;
t140 = t135+t137+t139;
t141 = 1.0./sqrt(t140);
t142 = t89.*5.65685424949238e2;
t143 = t90.^2;
t144 = t143+1.0e-8;
t145 = sqrt(t144);
t146 = t145.*5.65685424949238e2;
t252 = RL2.*5.65685424949238e2;
t147 = t142+t146-t252;
t148 = t92.^2;
t149 = t148.*2.0;
t150 = t93.^2;
t151 = t150.*2.0;
t152 = t94.^2;
t153 = t152.*2.0;
t154 = t149+t151+t153;
t155 = 1.0./sqrt(t154);
t156 = t101.*5.65685424949238e2;
t157 = t102.^2;
t158 = t157+1.0e-8;
t159 = sqrt(t158);
t160 = t159.*5.65685424949238e2;
t253 = RL1.*5.65685424949238e2;
t161 = t156+t160-t253;
t162 = t104.^2;
t163 = t162.*2.0;
t164 = t105.^2;
t165 = t164.*2.0;
t166 = t106.^2;
t167 = t166.*2.0;
t168 = t163+t165+t167;
t169 = 1.0./sqrt(t168);
t170 = t113.*5.65685424949238e2;
t171 = t114.^2;
t172 = t171+1.0e-8;
t173 = sqrt(t172);
t174 = t173.*5.65685424949238e2;
t254 = RL3.*5.65685424949238e2;
t175 = t170+t174-t254;
t176 = t116.^2;
t177 = t176.*2.0;
t178 = t117.^2;
t179 = t178.*2.0;
t180 = t118.^2;
t181 = t180.*2.0;
t182 = t177+t179+t181;
t183 = 1.0./sqrt(t182);
t184 = t125.*5.65685424949238e2;
t185 = t126.^2;
t186 = t185+1.0e-8;
t187 = sqrt(t186);
t188 = t187.*5.65685424949238e2;
t255 = RL4.*5.65685424949238e2;
t189 = t184+t188-t255;
t190 = t2+1.0e-4;
t191 = 1.0./sqrt(t190);
t192 = (p3.*t191)./2.0;
t193 = t192-1.0./2.0;
t194 = t4.^2;
t195 = t194.*2.0;
t196 = t5.^2;
t197 = t196.*2.0;
t198 = t6.^2;
t199 = t198.*2.0;
t200 = t195+t197+t199;
t201 = 1.0./sqrt(t200);
t202 = t13.*5.65685424949238e2;
t203 = t14.^2;
t204 = t203+1.0e-8;
t205 = sqrt(t204);
t206 = t205.*5.65685424949238e2;
t258 = RL5.*5.65685424949238e2;
t207 = t202+t206-t258;
t208 = t16.^2;
t209 = t208.*2.0;
t210 = t17.^2;
t211 = t210.*2.0;
t212 = t18.^2;
t213 = t212.*2.0;
t214 = t209+t211+t213;
t215 = 1.0./sqrt(t214);
t216 = t25.*5.65685424949238e2;
t217 = t26.^2;
t218 = t217+1.0e-8;
t219 = sqrt(t218);
t220 = t219.*5.65685424949238e2;
t259 = RL6.*5.65685424949238e2;
t221 = t216+t220-t259;
t222 = t28.^2;
t223 = t222.*2.0;
t224 = t29.^2;
t225 = t224.*2.0;
t226 = t30.^2;
t227 = t226.*2.0;
t228 = t223+t225+t227;
t229 = 1.0./sqrt(t228);
t230 = t37.*5.65685424949238e2;
t231 = t38.^2;
t232 = t231+1.0e-8;
t233 = sqrt(t232);
t234 = t233.*5.65685424949238e2;
t260 = RL7.*5.65685424949238e2;
t235 = t230+t234-t260;
t236 = t40.^2;
t237 = t236.*2.0;
t238 = t41.^2;
t239 = t238.*2.0;
t240 = t42.^2;
t241 = t240.*2.0;
t242 = t237+t239+t241;
t243 = 1.0./sqrt(t242);
t244 = t49.*5.65685424949238e2;
t245 = t50.^2;
t246 = t245+1.0e-8;
t247 = sqrt(t246);
t248 = t247.*5.65685424949238e2;
t261 = RL8.*5.65685424949238e2;
t249 = t244+t248-t261;
t250 = t56.*2.5e9;
t251 = t250+1.0e-4;
t256 = t2.*2.5e9;
t257 = t256+1.0e-4;
t262 = t7.*t201.*t207;
t263 = t19.*t215.*t221;
t264 = t31.*t229.*t235;
t265 = t43.*t243.*t249;
t328 = pDOT2.*t193.*5.0e1;
t266 = pDOT2+t262+t263+t264+t265-t328;
t267 = L1.^2;
t268 = (pDOT1.*t61)./1.0e4;
t269 = (pDOT2.*t64)./1.0e4;
t270 = (pDOT3.*t67)./1.0e4;
t271 = pDOT4.*t70;
t272 = pDOT5.*t73;
t273 = pDOT6.*t76;
t274 = t51.^2;
t275 = t274./1.0e4;
t276 = t77.^2;
t277 = t276./1.0e4;
t278 = t78.^2;
t279 = t278./1.0e4;
t318 = t267./1.0e4;
t319 = (pDOT4.*t61)./1.0e4;
t320 = (pDOT5.*t64)./1.0e4;
t321 = (pDOT6.*t67)./1.0e4;
t322 = pDOT1.*t70;
t323 = pDOT2.*t73;
t324 = pDOT3.*t76;
t280 = t268+t269+t270+t271+t272+t273+t275+t277+t279-t318-t319-t320-t321-t322-t323-t324;
t281 = t79.*t141.*t147;
t282 = t91.*t155.*t161;
t283 = t103.*t169.*t175;
t284 = t115.*t183.*t189;
t325 = pDOT4.*t133.*5.0e1;
t285 = pDOT4+t281+t282+t283+t284-t325;
t286 = t83.*t141.*t147;
t287 = t95.*t155.*t161;
t288 = t107.*t169.*t175;
t289 = t119.*t183.*t189;
t326 = pDOT5.*t133.*5.0e1;
t290 = pDOT5+t286+t287+t288+t289-t326;
t291 = t3.*t201.*t207;
t292 = t15.*t215.*t221;
t293 = t27.*t229.*t235;
t294 = t39.*t243.*t249;
t327 = pDOT1.*t193.*5.0e1;
t295 = pDOT1+t291+t292+t293+t294-t327;
t296 = p6.*2.5e4;
t297 = 1.0./sqrt(t251);
t298 = p6.*t297.*2.5e4;
t299 = t298-1.0./2.0;
t300 = sqrt(t251);
t301 = t84.*t141.*t147;
t302 = t96.*t155.*t161;
t303 = t108.*t169.*t175;
t304 = t120.*t183.*t189;
t329 = pDOT6.*t299.*3.0e2;
t330 = t300./2.0;
t305 = pDOT6+t296+t301+t302+t303+t304-t329-t330+1.4715;
t306 = p3.*2.5e4;
t307 = 1.0./sqrt(t257);
t308 = p3.*t307.*2.5e4;
t309 = t308-1.0./2.0;
t310 = sqrt(t257);
t311 = t8.*t201.*t207;
t312 = t20.*t215.*t221;
t313 = t32.*t229.*t235;
t314 = t44.*t243.*t249;
t316 = pDOT3.*t309.*3.0e2;
t317 = t310./2.0;
t315 = pDOT3+t306+t311+t312+t313+t314-t316-t317+1.4715;
t331 = (t51.*t58.*t61)./4.0;
t332 = t331-1.0;
t333 = (t51.*t58.*t280)./4.0;
t334 = t51.*t58.*t64.*t266.*(5.0./3.0);
t335 = t51.*t58.*t67.*t315.*(5.0./3.0);
t336 = (t58.*t64.*t77)./4.0;
t337 = t336-1.0;
t338 = (t58.*t77.*t280)./4.0;
t339 = t58.*t61.*t77.*t295.*(5.0./3.0);
t340 = t58.*t67.*t77.*t315.*(5.0./3.0);
t341 = (t58.*t67.*t78)./4.0;
t342 = t341-1.0;
t343 = (t58.*t78.*t280)./4.0;
t344 = t58.*t61.*t78.*t295.*(5.0./3.0);
t345 = t58.*t64.*t78.*t266.*(5.0./3.0);
t346 = p9.^2;
t347 = p7-p19;
t348 = abs(t347);
t351 = p8-p20;
t349 = abs(t351);
t352 = p9-p21;
t350 = abs(t352);
t353 = t347.^2;
t354 = t351.^2;
t355 = t352.^2;
t356 = t353+t354+t355;
t357 = sqrt(t356);
t358 = RL9-t357;
t359 = p7-p22;
t360 = abs(t359);
t363 = p8-p23;
t361 = abs(t363);
t364 = p9-p24;
t362 = abs(t364);
t365 = t359.^2;
t366 = t363.^2;
t367 = t364.^2;
t368 = t365+t366+t367;
t369 = sqrt(t368);
t370 = RL10-t369;
t371 = p7-p25;
t372 = abs(t371);
t375 = p8-p26;
t373 = abs(t375);
t376 = p9-p27;
t374 = abs(t376);
t377 = t371.^2;
t378 = t375.^2;
t379 = t376.^2;
t380 = t377+t378+t379;
t381 = sqrt(t380);
t382 = RL12-t381;
t383 = p7-p31;
t384 = abs(t383);
t387 = p8-p32;
t385 = abs(t387);
t388 = p9-p33;
t386 = abs(t388);
t389 = t383.^2;
t390 = t387.^2;
t391 = t388.^2;
t392 = t389+t390+t391;
t393 = sqrt(t392);
t394 = RL11-t393;
t395 = p7-p10;
t396 = p7.^2;
t397 = p8.^2;
t398 = p10.^2;
t399 = p11.^2;
t400 = p12.^2;
t483 = p7.*p10.*2.0;
t484 = p8.*p11.*2.0;
t485 = p9.*p12.*2.0;
t401 = t346+t396+t397+t398+t399+t400-t483-t484-t485;
t402 = 1.0./t401;
t403 = p7.*2.0;
t404 = p10.*2.0;
t405 = t403-t404;
t406 = p8.*2.0;
t407 = p11.*2.0;
t408 = t406-t407;
t409 = p9.*2.0;
t410 = p12.*2.0;
t411 = t409-t410;
t412 = pDOT7.*2.0;
t413 = pDOT10.*2.0;
t414 = t412-t413;
t415 = pDOT8.*2.0;
t416 = pDOT11.*2.0;
t417 = t415-t416;
t418 = pDOT9.*2.0;
t419 = pDOT12.*2.0;
t420 = t418-t419;
t421 = p8-p11;
t422 = p9-p12;
t423 = t346+1.0e-4;
t424 = 1.0./sqrt(t423);
t425 = (p9.*t424)./2.0;
t426 = t425-1.0./2.0;
t427 = t348.^2;
t428 = t427.*2.0;
t429 = t349.^2;
t430 = t429.*2.0;
t431 = t350.^2;
t432 = t431.*2.0;
t433 = t428+t430+t432;
t434 = 1.0./sqrt(t433);
t435 = t357.*5.65685424949238e2;
t436 = t358.^2;
t437 = t436+1.0e-8;
t438 = sqrt(t437);
t439 = t438.*5.65685424949238e2;
t596 = RL9.*5.65685424949238e2;
t440 = t435+t439-t596;
t441 = t360.^2;
t442 = t441.*2.0;
t443 = t361.^2;
t444 = t443.*2.0;
t445 = t362.^2;
t446 = t445.*2.0;
t447 = t442+t444+t446;
t448 = 1.0./sqrt(t447);
t449 = t369.*5.65685424949238e2;
t450 = t370.^2;
t451 = t450+1.0e-8;
t452 = sqrt(t451);
t453 = t452.*5.65685424949238e2;
t597 = RL10.*5.65685424949238e2;
t454 = t449+t453-t597;
t455 = t372.^2;
t456 = t455.*2.0;
t457 = t373.^2;
t458 = t457.*2.0;
t459 = t374.^2;
t460 = t459.*2.0;
t461 = t456+t458+t460;
t462 = 1.0./sqrt(t461);
t463 = t381.*5.65685424949238e2;
t464 = t382.^2;
t465 = t464+1.0e-8;
t466 = sqrt(t465);
t467 = t466.*5.65685424949238e2;
t598 = RL12.*5.65685424949238e2;
t468 = t463+t467-t598;
t469 = t384.^2;
t470 = t469.*2.0;
t471 = t385.^2;
t472 = t471.*2.0;
t473 = t386.^2;
t474 = t473.*2.0;
t475 = t470+t472+t474;
t476 = 1.0./sqrt(t475);
t477 = t393.*5.65685424949238e2;
t478 = t394.^2;
t479 = t478+1.0e-8;
t480 = sqrt(t479);
t481 = t480.*5.65685424949238e2;
t599 = RL11.*5.65685424949238e2;
t482 = t477+t481-t599;
t486 = p10-p13;
t487 = abs(t486);
t490 = p11-p14;
t488 = abs(t490);
t491 = p12-p15;
t489 = abs(t491);
t492 = t486.^2;
t493 = t490.^2;
t494 = t491.^2;
t495 = t492+t493+t494;
t496 = sqrt(t495);
t497 = RL14-t496;
t498 = p10-p16;
t499 = abs(t498);
t502 = p11-p17;
t500 = abs(t502);
t503 = p12-p18;
t501 = abs(t503);
t504 = t498.^2;
t505 = t502.^2;
t506 = t503.^2;
t507 = t504+t505+t506;
t508 = sqrt(t507);
t509 = RL13-t508;
t510 = p10-p25;
t511 = abs(t510);
t514 = p11-p26;
t512 = abs(t514);
t515 = p12-p27;
t513 = abs(t515);
t516 = t510.^2;
t517 = t514.^2;
t518 = t515.^2;
t519 = t516+t517+t518;
t520 = sqrt(t519);
t521 = RL16-t520;
t522 = p10-p31;
t523 = abs(t522);
t526 = p11-p32;
t524 = abs(t526);
t527 = p12-p33;
t525 = abs(t527);
t528 = t522.^2;
t529 = t526.^2;
t530 = t527.^2;
t531 = t528+t529+t530;
t532 = sqrt(t531);
t533 = RL15-t532;
t534 = t400+1.0e-4;
t535 = 1.0./sqrt(t534);
t536 = (p12.*t535)./2.0;
t537 = t536-1.0./2.0;
t538 = t487.^2;
t539 = t538.*2.0;
t540 = t488.^2;
t541 = t540.*2.0;
t542 = t489.^2;
t543 = t542.*2.0;
t544 = t539+t541+t543;
t545 = 1.0./sqrt(t544);
t546 = t496.*5.65685424949238e2;
t547 = t497.^2;
t548 = t547+1.0e-8;
t549 = sqrt(t548);
t550 = t549.*5.65685424949238e2;
t602 = RL14.*5.65685424949238e2;
t551 = t546+t550-t602;
t552 = t499.^2;
t553 = t552.*2.0;
t554 = t500.^2;
t555 = t554.*2.0;
t556 = t501.^2;
t557 = t556.*2.0;
t558 = t553+t555+t557;
t559 = 1.0./sqrt(t558);
t560 = t508.*5.65685424949238e2;
t561 = t509.^2;
t562 = t561+1.0e-8;
t563 = sqrt(t562);
t564 = t563.*5.65685424949238e2;
t603 = RL13.*5.65685424949238e2;
t565 = t560+t564-t603;
t566 = t511.^2;
t567 = t566.*2.0;
t568 = t512.^2;
t569 = t568.*2.0;
t570 = t513.^2;
t571 = t570.*2.0;
t572 = t567+t569+t571;
t573 = 1.0./sqrt(t572);
t574 = t520.*5.65685424949238e2;
t575 = t521.^2;
t576 = t575+1.0e-8;
t577 = sqrt(t576);
t578 = t577.*5.65685424949238e2;
t604 = RL16.*5.65685424949238e2;
t579 = t574+t578-t604;
t580 = t523.^2;
t581 = t580.*2.0;
t582 = t524.^2;
t583 = t582.*2.0;
t584 = t525.^2;
t585 = t584.*2.0;
t586 = t581+t583+t585;
t587 = 1.0./sqrt(t586);
t588 = t532.*5.65685424949238e2;
t589 = t533.^2;
t590 = t589+1.0e-8;
t591 = sqrt(t590);
t592 = t591.*5.65685424949238e2;
t605 = RL15.*5.65685424949238e2;
t593 = t588+t592-t605;
t594 = t346.*2.5e9;
t595 = t594+1.0e-4;
t600 = t400.*2.5e9;
t601 = t600+1.0e-4;
t606 = t351.*t434.*t440;
t607 = t363.*t448.*t454;
t608 = t375.*t462.*t468;
t609 = t387.*t476.*t482;
t670 = pDOT8.*t426.*5.0e1;
t610 = pDOT8+t606+t607+t608+t609-t670;
t611 = L2.^2;
t612 = (pDOT7.*t405)./1.0e4;
t613 = (pDOT8.*t408)./1.0e4;
t614 = (pDOT9.*t411)./1.0e4;
t615 = pDOT10.*t414;
t616 = pDOT11.*t417;
t617 = pDOT12.*t420;
t618 = t395.^2;
t619 = t618./1.0e4;
t620 = t421.^2;
t621 = t620./1.0e4;
t622 = t422.^2;
t623 = t622./1.0e4;
t662 = t611./1.0e4;
t663 = (pDOT10.*t405)./1.0e4;
t664 = (pDOT11.*t408)./1.0e4;
t665 = (pDOT12.*t411)./1.0e4;
t666 = pDOT7.*t414;
t667 = pDOT8.*t417;
t668 = pDOT9.*t420;
t624 = t612+t613+t614+t615+t616+t617+t619+t621+t623-t662-t663-t664-t665-t666-t667-t668;
t625 = t347.*t434.*t440;
t626 = t359.*t448.*t454;
t627 = t371.*t462.*t468;
t628 = t383.*t476.*t482;
t669 = pDOT7.*t426.*5.0e1;
t629 = pDOT7+t625+t626+t627+t628-t669;
t630 = t486.*t545.*t551;
t631 = t498.*t559.*t565;
t632 = t510.*t573.*t579;
t633 = t522.*t587.*t593;
t671 = pDOT10.*t537.*5.0e1;
t634 = pDOT10+t630+t631+t632+t633-t671;
t635 = t490.*t545.*t551;
t636 = t502.*t559.*t565;
t637 = t514.*t573.*t579;
t638 = t526.*t587.*t593;
t672 = pDOT11.*t537.*5.0e1;
t639 = pDOT11+t635+t636+t637+t638-t672;
t640 = p9.*2.5e4;
t641 = 1.0./sqrt(t595);
t642 = p9.*t641.*2.5e4;
t643 = t642-1.0./2.0;
t644 = sqrt(t595);
t645 = t352.*t434.*t440;
t646 = t364.*t448.*t454;
t647 = t376.*t462.*t468;
t648 = t388.*t476.*t482;
t660 = pDOT9.*t643.*3.0e2;
t661 = t644./2.0;
t649 = pDOT9+t640+t645+t646+t647+t648-t660-t661+1.4715;
t650 = p12.*2.5e4;
t651 = 1.0./sqrt(t601);
t652 = p12.*t651.*2.5e4;
t653 = t652-1.0./2.0;
t654 = sqrt(t601);
t655 = t491.*t545.*t551;
t656 = t503.*t559.*t565;
t657 = t515.*t573.*t579;
t658 = t527.*t587.*t593;
t673 = pDOT12.*t653.*3.0e2;
t674 = t654./2.0;
t659 = pDOT12+t650+t655+t656+t657+t658-t673-t674+1.4715;
t675 = (t395.*t402.*t405)./4.0;
t676 = t675-1.0;
t677 = (t395.*t402.*t624)./4.0;
t678 = t395.*t402.*t408.*t610.*(5.0./3.0);
t679 = t395.*t402.*t411.*t649.*(5.0./3.0);
t680 = (t402.*t408.*t421)./4.0;
t681 = t680-1.0;
t682 = (t402.*t421.*t624)./4.0;
t683 = t402.*t405.*t421.*t629.*(5.0./3.0);
t684 = t402.*t411.*t421.*t649.*(5.0./3.0);
t685 = (t402.*t411.*t422)./4.0;
t686 = t685-1.0;
t687 = (t402.*t422.*t624)./4.0;
t688 = t402.*t405.*t422.*t629.*(5.0./3.0);
t689 = t402.*t408.*t422.*t610.*(5.0./3.0);
t690 = p15.^2;
t691 = p13-p31;
t692 = abs(t691);
t695 = p14-p32;
t693 = abs(t695);
t696 = p15-p33;
t694 = abs(t696);
t697 = t691.^2;
t698 = t695.^2;
t699 = t696.^2;
t700 = t697+t698+t699;
t701 = sqrt(t700);
t702 = RL20-t701;
t703 = p13-p34;
t704 = abs(t703);
t707 = p14-p35;
t705 = abs(t707);
t708 = p15-p36;
t706 = abs(t708);
t709 = t703.^2;
t710 = t707.^2;
t711 = t708.^2;
t712 = t709+t710+t711;
t713 = sqrt(t712);
t714 = RL19-t713;
t715 = p13-p16;
t716 = p13.^2;
t717 = p14.^2;
t718 = p16.^2;
t719 = p17.^2;
t720 = p18.^2;
t767 = p13.*p16.*2.0;
t768 = p14.*p17.*2.0;
t769 = p15.*p18.*2.0;
t721 = t690+t716+t717+t718+t719+t720-t767-t768-t769;
t722 = 1.0./t721;
t723 = p13.*2.0;
t724 = p16.*2.0;
t725 = t723-t724;
t726 = p14.*2.0;
t727 = p17.*2.0;
t728 = t726-t727;
t729 = p15.*2.0;
t730 = p18.*2.0;
t731 = t729-t730;
t732 = pDOT13.*2.0;
t733 = pDOT16.*2.0;
t734 = t732-t733;
t735 = pDOT14.*2.0;
t736 = pDOT17.*2.0;
t737 = t735-t736;
t738 = pDOT15.*2.0;
t739 = pDOT18.*2.0;
t740 = t738-t739;
t741 = p14-p17;
t742 = p15-p18;
t743 = p16-p25;
t744 = abs(t743);
t747 = p17-p26;
t745 = abs(t747);
t748 = p18-p27;
t746 = abs(t748);
t749 = t743.^2;
t750 = t747.^2;
t751 = t748.^2;
t752 = t749+t750+t751;
t753 = sqrt(t752);
t754 = RL18-t753;
t755 = p16-p28;
t756 = abs(t755);
t759 = p17-p29;
t757 = abs(t759);
t760 = p18-p30;
t758 = abs(t760);
t761 = t755.^2;
t762 = t759.^2;
t763 = t760.^2;
t764 = t761+t762+t763;
t765 = sqrt(t764);
t766 = RL17-t765;
t770 = t690+1.0e-4;
t771 = 1.0./sqrt(t770);
t772 = (p15.*t771)./2.0;
t773 = t772-1.0./2.0;
t774 = t692.^2;
t775 = t774.*2.0;
t776 = t693.^2;
t777 = t776.*2.0;
t778 = t694.^2;
t779 = t778.*2.0;
t780 = t775+t777+t779;
t781 = 1.0./sqrt(t780);
t782 = t701.*5.65685424949238e2;
t783 = t702.^2;
t784 = t783+1.0e-8;
t785 = sqrt(t784);
t786 = t785.*5.65685424949238e2;
t836 = RL20.*5.65685424949238e2;
t787 = t782+t786-t836;
t788 = t704.^2;
t789 = t788.*2.0;
t790 = t705.^2;
t791 = t790.*2.0;
t792 = t706.^2;
t793 = t792.*2.0;
t794 = t789+t791+t793;
t795 = 1.0./sqrt(t794);
t796 = t713.*5.65685424949238e2;
t797 = t714.^2;
t798 = t797+1.0e-8;
t799 = sqrt(t798);
t800 = t799.*5.65685424949238e2;
t837 = RL19.*5.65685424949238e2;
t801 = t796+t800-t837;
t802 = t720+1.0e-4;
t803 = 1.0./sqrt(t802);
t804 = (p18.*t803)./2.0;
t805 = t804-1.0./2.0;
t806 = t744.^2;
t807 = t806.*2.0;
t808 = t745.^2;
t809 = t808.*2.0;
t810 = t746.^2;
t811 = t810.*2.0;
t812 = t807+t809+t811;
t813 = 1.0./sqrt(t812);
t814 = t753.*5.65685424949238e2;
t815 = t754.^2;
t816 = t815+1.0e-8;
t817 = sqrt(t816);
t818 = t817.*5.65685424949238e2;
t840 = RL18.*5.65685424949238e2;
t819 = t814+t818-t840;
t820 = t756.^2;
t821 = t820.*2.0;
t822 = t757.^2;
t823 = t822.*2.0;
t824 = t758.^2;
t825 = t824.*2.0;
t826 = t821+t823+t825;
t827 = 1.0./sqrt(t826);
t828 = t765.*5.65685424949238e2;
t829 = t766.^2;
t830 = t829+1.0e-8;
t831 = sqrt(t830);
t832 = t831.*5.65685424949238e2;
t841 = RL17.*5.65685424949238e2;
t833 = t828+t832-t841;
t834 = t690.*2.5e9;
t835 = t834+1.0e-4;
t838 = t720.*2.5e9;
t839 = t838+1.0e-4;
t842 = t695.*t781.*t787;
t843 = t707.*t795.*t801;
t895 = pDOT14.*t773.*5.0e1;
t844 = pDOT14-t286-t635+t842+t843-t895;
t845 = L3.^2;
t846 = (pDOT13.*t725)./1.0e4;
t847 = (pDOT14.*t728)./1.0e4;
t848 = (pDOT15.*t731)./1.0e4;
t849 = pDOT16.*t734;
t850 = pDOT17.*t737;
t851 = pDOT18.*t740;
t852 = t715.^2;
t853 = t852./1.0e4;
t854 = t741.^2;
t855 = t854./1.0e4;
t856 = t742.^2;
t857 = t856./1.0e4;
t886 = t845./1.0e4;
t887 = (pDOT16.*t725)./1.0e4;
t888 = (pDOT17.*t728)./1.0e4;
t889 = (pDOT18.*t731)./1.0e4;
t890 = pDOT13.*t734;
t891 = pDOT14.*t737;
t892 = pDOT15.*t740;
t858 = t846+t847+t848+t849+t850+t851+t853+t855+t857-t886-t887-t888-t889-t890-t891-t892;
t859 = t691.*t781.*t787;
t860 = t703.*t795.*t801;
t893 = pDOT13.*t773.*5.0e1;
t861 = pDOT13-t281-t630+t859+t860-t893;
t862 = t743.*t813.*t819;
t863 = t755.*t827.*t833;
t894 = pDOT16.*t805.*5.0e1;
t864 = pDOT16-t282-t631+t862+t863-t894;
t865 = t747.*t813.*t819;
t866 = t759.*t827.*t833;
t896 = pDOT17.*t805.*5.0e1;
t867 = pDOT17-t287-t636+t865+t866-t896;
t868 = p15.*2.5e4;
t869 = 1.0./sqrt(t835);
t870 = p15.*t869.*2.5e4;
t871 = t870-1.0./2.0;
t872 = sqrt(t835);
t873 = t696.*t781.*t787;
t874 = t708.*t795.*t801;
t884 = pDOT15.*t871.*3.0e2;
t885 = t872./2.0;
t875 = pDOT15-t301-t655+t868+t873+t874-t884-t885+1.4715;
t876 = p18.*2.5e4;
t877 = 1.0./sqrt(t839);
t878 = p18.*t877.*2.5e4;
t879 = t878-1.0./2.0;
t880 = sqrt(t839);
t881 = t748.*t813.*t819;
t882 = t760.*t827.*t833;
t897 = pDOT18.*t879.*3.0e2;
t898 = t880./2.0;
t883 = pDOT18-t302-t656+t876+t881+t882-t897-t898+1.4715;
t899 = (t715.*t722.*t725)./4.0;
t900 = t899-1.0;
t901 = (t715.*t722.*t858)./4.0;
t902 = t715.*t722.*t728.*t844.*(5.0./3.0);
t903 = t715.*t722.*t731.*t875.*(5.0./3.0);
t904 = (t722.*t728.*t741)./4.0;
t905 = t904-1.0;
t906 = (t722.*t741.*t858)./4.0;
t907 = t722.*t725.*t741.*t861.*(5.0./3.0);
t908 = t722.*t731.*t741.*t875.*(5.0./3.0);
t909 = (t722.*t731.*t742)./4.0;
t910 = t909-1.0;
t911 = (t722.*t742.*t858)./4.0;
t912 = t722.*t725.*t742.*t861.*(5.0./3.0);
t913 = t722.*t728.*t742.*t844.*(5.0./3.0);
t914 = p21.^2;
t915 = p19-p31;
t916 = abs(t915);
t919 = p20-p32;
t917 = abs(t919);
t920 = p21-p33;
t918 = abs(t920);
t921 = t915.^2;
t922 = t919.^2;
t923 = t920.^2;
t924 = t921+t922+t923;
t925 = sqrt(t924);
t926 = RL22-t925;
t927 = p19-p34;
t928 = abs(t927);
t931 = p20-p35;
t929 = abs(t931);
t932 = p21-p36;
t930 = abs(t932);
t933 = t927.^2;
t934 = t931.^2;
t935 = t932.^2;
t936 = t933+t934+t935;
t937 = sqrt(t936);
t938 = RL21-t937;
t939 = p19-p22;
t940 = p19.^2;
t941 = p20.^2;
t942 = p22.^2;
t943 = p23.^2;
t944 = p24.^2;
t999 = p19.*p22.*2.0;
t1000 = p20.*p23.*2.0;
t1001 = p21.*p24.*2.0;
t945 = t914+t940+t941+t942+t943+t944-t999-t1000-t1001;
t946 = 1.0./t945;
t947 = p19.*2.0;
t948 = p22.*2.0;
t949 = t947-t948;
t950 = p20.*2.0;
t951 = p23.*2.0;
t952 = t950-t951;
t953 = p21.*2.0;
t954 = p24.*2.0;
t955 = t953-t954;
t956 = pDOT19.*2.0;
t957 = pDOT22.*2.0;
t958 = t956-t957;
t959 = pDOT20.*2.0;
t960 = pDOT23.*2.0;
t961 = t959-t960;
t962 = pDOT21.*2.0;
t963 = pDOT24.*2.0;
t964 = t962-t963;
t965 = p20-p23;
t966 = p21-p24;
t967 = t914+1.0e-4;
t968 = 1.0./sqrt(t967);
t969 = (p21.*t968)./2.0;
t970 = t969-1.0./2.0;
t971 = t916.^2;
t972 = t971.*2.0;
t973 = t917.^2;
t974 = t973.*2.0;
t975 = t918.^2;
t976 = t975.*2.0;
t977 = t972+t974+t976;
t978 = 1.0./sqrt(t977);
t979 = t925.*5.65685424949238e2;
t980 = t926.^2;
t981 = t980+1.0e-8;
t982 = sqrt(t981);
t983 = t982.*5.65685424949238e2;
t1060 = RL22.*5.65685424949238e2;
t984 = t979+t983-t1060;
t985 = t928.^2;
t986 = t985.*2.0;
t987 = t929.^2;
t988 = t987.*2.0;
t989 = t930.^2;
t990 = t989.*2.0;
t991 = t986+t988+t990;
t992 = 1.0./sqrt(t991);
t993 = t937.*5.65685424949238e2;
t994 = t938.^2;
t995 = t994+1.0e-8;
t996 = sqrt(t995);
t997 = t996.*5.65685424949238e2;
t1061 = RL21.*5.65685424949238e2;
t998 = t993+t997-t1061;
t1002 = p22-p25;
t1003 = abs(t1002);
t1006 = p23-p26;
t1004 = abs(t1006);
t1007 = p24-p27;
t1005 = abs(t1007);
t1008 = t1002.^2;
t1009 = t1006.^2;
t1010 = t1007.^2;
t1011 = t1008+t1009+t1010;
t1012 = sqrt(t1011);
t1013 = RL24-t1012;
t1014 = p22-p28;
t1015 = abs(t1014);
t1018 = p23-p29;
t1016 = abs(t1018);
t1019 = p24-p30;
t1017 = abs(t1019);
t1020 = t1014.^2;
t1021 = t1018.^2;
t1022 = t1019.^2;
t1023 = t1020+t1021+t1022;
t1024 = sqrt(t1023);
t1025 = RL23-t1024;
t1026 = t944+1.0e-4;
t1027 = 1.0./sqrt(t1026);
t1028 = (p24.*t1027)./2.0;
t1029 = t1028-1.0./2.0;
t1030 = t1003.^2;
t1031 = t1030.*2.0;
t1032 = t1004.^2;
t1033 = t1032.*2.0;
t1034 = t1005.^2;
t1035 = t1034.*2.0;
t1036 = t1031+t1033+t1035;
t1037 = 1.0./sqrt(t1036);
t1038 = t1012.*5.65685424949238e2;
t1039 = t1013.^2;
t1040 = t1039+1.0e-8;
t1041 = sqrt(t1040);
t1042 = t1041.*5.65685424949238e2;
t1064 = RL24.*5.65685424949238e2;
t1043 = t1038+t1042-t1064;
t1044 = t1015.^2;
t1045 = t1044.*2.0;
t1046 = t1016.^2;
t1047 = t1046.*2.0;
t1048 = t1017.^2;
t1049 = t1048.*2.0;
t1050 = t1045+t1047+t1049;
t1051 = 1.0./sqrt(t1050);
t1052 = t1024.*5.65685424949238e2;
t1053 = t1025.^2;
t1054 = t1053+1.0e-8;
t1055 = sqrt(t1054);
t1056 = t1055.*5.65685424949238e2;
t1065 = RL23.*5.65685424949238e2;
t1057 = t1052+t1056-t1065;
t1058 = t914.*2.5e9;
t1059 = t1058+1.0e-4;
t1062 = t944.*2.5e9;
t1063 = t1062+1.0e-4;
t1066 = t919.*t978.*t984;
t1067 = t931.*t992.*t998;
t1118 = pDOT20.*t970.*5.0e1;
t1068 = pDOT20-t262-t606+t1066+t1067-t1118;
t1069 = L4.^2;
t1070 = (pDOT19.*t949)./1.0e4;
t1071 = (pDOT20.*t952)./1.0e4;
t1072 = (pDOT21.*t955)./1.0e4;
t1073 = pDOT22.*t958;
t1074 = pDOT23.*t961;
t1075 = pDOT24.*t964;
t1076 = t939.^2;
t1077 = t1076./1.0e4;
t1078 = t965.^2;
t1079 = t1078./1.0e4;
t1080 = t966.^2;
t1081 = t1080./1.0e4;
t1110 = t1069./1.0e4;
t1111 = (pDOT22.*t949)./1.0e4;
t1112 = (pDOT23.*t952)./1.0e4;
t1113 = (pDOT24.*t955)./1.0e4;
t1114 = pDOT19.*t958;
t1115 = pDOT20.*t961;
t1116 = pDOT21.*t964;
t1082 = t1070+t1071+t1072+t1073+t1074+t1075+t1077+t1079+t1081-t1110-t1111-t1112-t1113-t1114-t1115-t1116;
t1083 = t915.*t978.*t984;
t1084 = t927.*t992.*t998;
t1117 = pDOT19.*t970.*5.0e1;
t1085 = pDOT19-t291-t625+t1083+t1084-t1117;
t1086 = t1002.*t1037.*t1043;
t1087 = t1014.*t1051.*t1057;
t1119 = pDOT22.*t1029.*5.0e1;
t1088 = pDOT22-t292-t626+t1086+t1087-t1119;
t1089 = t1006.*t1037.*t1043;
t1090 = t1018.*t1051.*t1057;
t1120 = pDOT23.*t1029.*5.0e1;
t1091 = pDOT23-t263-t607+t1089+t1090-t1120;
t1092 = p21.*2.5e4;
t1093 = 1.0./sqrt(t1059);
t1094 = p21.*t1093.*2.5e4;
t1095 = t1094-1.0./2.0;
t1096 = sqrt(t1059);
t1097 = t920.*t978.*t984;
t1098 = t932.*t992.*t998;
t1108 = pDOT21.*t1095.*3.0e2;
t1109 = t1096./2.0;
t1099 = pDOT21-t311-t645+t1092+t1097+t1098-t1108-t1109+1.4715;
t1100 = p24.*2.5e4;
t1101 = 1.0./sqrt(t1063);
t1102 = p24.*t1101.*2.5e4;
t1103 = t1102-1.0./2.0;
t1104 = sqrt(t1063);
t1105 = t1007.*t1037.*t1043;
t1106 = t1019.*t1051.*t1057;
t1121 = pDOT24.*t1103.*3.0e2;
t1122 = t1104./2.0;
t1107 = pDOT24-t312-t646+t1100+t1105+t1106-t1121-t1122+1.4715;
t1123 = (t939.*t946.*t949)./4.0;
t1124 = t1123-1.0;
t1125 = (t939.*t946.*t1082)./4.0;
t1126 = t939.*t946.*t952.*t1068.*(5.0./3.0);
t1127 = t939.*t946.*t955.*t1099.*(5.0./3.0);
t1128 = (t946.*t952.*t965)./4.0;
t1129 = t1128-1.0;
t1130 = (t946.*t965.*t1082)./4.0;
t1131 = t946.*t949.*t965.*t1085.*(5.0./3.0);
t1132 = t946.*t955.*t965.*t1099.*(5.0./3.0);
t1133 = (t946.*t955.*t966)./4.0;
t1134 = t1133-1.0;
t1135 = (t946.*t966.*t1082)./4.0;
t1136 = t946.*t949.*t966.*t1085.*(5.0./3.0);
t1137 = t946.*t952.*t966.*t1068.*(5.0./3.0);
t1138 = p27.^2;
t1139 = p25-p28;
t1140 = p25.^2;
t1141 = p26.^2;
t1142 = p28.^2;
t1143 = p29.^2;
t1144 = p30.^2;
t1169 = p25.*p28.*2.0;
t1170 = p26.*p29.*2.0;
t1171 = p27.*p30.*2.0;
t1145 = t1138+t1140+t1141+t1142+t1143+t1144-t1169-t1170-t1171;
t1146 = 1.0./t1145;
t1147 = p25.*2.0;
t1148 = p28.*2.0;
t1149 = t1147-t1148;
t1150 = p26.*2.0;
t1151 = p29.*2.0;
t1152 = t1150-t1151;
t1153 = p27.*2.0;
t1154 = p30.*2.0;
t1155 = t1153-t1154;
t1156 = pDOT25.*2.0;
t1157 = pDOT28.*2.0;
t1158 = t1156-t1157;
t1159 = pDOT26.*2.0;
t1160 = pDOT29.*2.0;
t1161 = t1159-t1160;
t1162 = pDOT27.*2.0;
t1163 = pDOT30.*2.0;
t1164 = t1162-t1163;
t1165 = p26-p29;
t1166 = p27-p30;
t1167 = t1144.*2.5e9;
t1168 = t1167+1.0e-4;
t1172 = t1138.*2.5e9;
t1173 = t1172+1.0e-4;
t1174 = t1144+1.0e-4;
t1175 = 1.0./sqrt(t1174);
t1176 = (p30.*t1175)./2.0;
t1177 = t1176-1.0./2.0;
t1178 = t1138+1.0e-4;
t1179 = 1.0./sqrt(t1178);
t1180 = (p27.*t1179)./2.0;
t1181 = t1180-1.0./2.0;
t1182 = pDOT26.*t1181.*5.0e1;
t1183 = -pDOT26+t608+t637+t865+t1089+t1182;
t1184 = L5.^2;
t1185 = (pDOT25.*t1149)./1.0e4;
t1186 = (pDOT26.*t1152)./1.0e4;
t1187 = (pDOT27.*t1155)./1.0e4;
t1188 = pDOT28.*t1158;
t1189 = pDOT29.*t1161;
t1190 = pDOT30.*t1164;
t1191 = t1139.^2;
t1192 = t1191./1.0e4;
t1193 = t1165.^2;
t1194 = t1193./1.0e4;
t1195 = t1166.^2;
t1196 = t1195./1.0e4;
t1219 = t1184./1.0e4;
t1220 = (pDOT28.*t1149)./1.0e4;
t1221 = (pDOT29.*t1152)./1.0e4;
t1222 = (pDOT30.*t1155)./1.0e4;
t1223 = pDOT25.*t1158;
t1224 = pDOT26.*t1161;
t1225 = pDOT27.*t1164;
t1197 = t1185+t1186+t1187+t1188+t1189+t1190+t1192+t1194+t1196-t1219-t1220-t1221-t1222-t1223-t1224-t1225;
t1198 = 1.0./sqrt(t1168);
t1199 = p30.*t1198.*2.5e4;
t1200 = t1199-1.0./2.0;
t1201 = pDOT30.*t1200.*3.0e2;
t1202 = sqrt(t1168);
t1203 = t1202./2.0;
t1226 = p30.*2.5e4;
t1204 = -pDOT30+t303+t313+t882+t1106+t1201+t1203-t1226-1.4715;
t1205 = 1.0./sqrt(t1173);
t1206 = p27.*t1205.*2.5e4;
t1207 = t1206-1.0./2.0;
t1208 = pDOT27.*t1207.*3.0e2;
t1209 = sqrt(t1173);
t1210 = t1209./2.0;
t1218 = p27.*2.5e4;
t1211 = -pDOT27+t647+t657+t881+t1105+t1208+t1210-t1218-1.4715;
t1212 = pDOT28.*t1177.*5.0e1;
t1213 = -pDOT28+t283+t293+t863+t1087+t1212;
t1214 = pDOT29.*t1177.*5.0e1;
t1215 = -pDOT29+t264+t288+t866+t1090+t1214;
t1216 = pDOT25.*t1181.*5.0e1;
t1217 = -pDOT25+t627+t632+t862+t1086+t1216;
t1227 = (t1139.*t1146.*t1149)./4.0;
t1228 = t1227-1.0;
t1229 = (t1139.*t1146.*t1197)./4.0;
t1230 = t1139.*t1146.*t1155.*t1204.*(5.0./3.0);
t1231 = t1139.*t1146.*t1152.*t1215.*(5.0./3.0);
t1232 = (t1146.*t1152.*t1165)./4.0;
t1233 = t1232-1.0;
t1234 = (t1146.*t1165.*t1197)./4.0;
t1235 = t1146.*t1155.*t1165.*t1204.*(5.0./3.0);
t1236 = t1146.*t1149.*t1165.*t1213.*(5.0./3.0);
t1237 = (t1146.*t1155.*t1166)./4.0;
t1238 = t1237-1.0;
t1239 = (t1146.*t1166.*t1197)./4.0;
t1240 = t1146.*t1149.*t1166.*t1213.*(5.0./3.0);
t1241 = t1146.*t1152.*t1166.*t1215.*(5.0./3.0);
t1242 = p33.^2;
t1243 = p31-p34;
t1244 = p31.^2;
t1245 = p32.^2;
t1246 = p34.^2;
t1247 = p35.^2;
t1248 = p36.^2;
t1273 = p31.*p34.*2.0;
t1274 = p32.*p35.*2.0;
t1275 = p33.*p36.*2.0;
t1249 = t1242+t1244+t1245+t1246+t1247+t1248-t1273-t1274-t1275;
t1250 = 1.0./t1249;
t1251 = p31.*2.0;
t1252 = p34.*2.0;
t1253 = t1251-t1252;
t1254 = p32.*2.0;
t1255 = p35.*2.0;
t1256 = t1254-t1255;
t1257 = p33.*2.0;
t1258 = p36.*2.0;
t1259 = t1257-t1258;
t1260 = pDOT31.*2.0;
t1261 = pDOT34.*2.0;
t1262 = t1260-t1261;
t1263 = pDOT32.*2.0;
t1264 = pDOT35.*2.0;
t1265 = t1263-t1264;
t1266 = pDOT33.*2.0;
t1267 = pDOT36.*2.0;
t1268 = t1266-t1267;
t1269 = p32-p35;
t1270 = p33-p36;
t1271 = t1248.*2.5e9;
t1272 = t1271+1.0e-4;
t1276 = t1242.*2.5e9;
t1277 = t1276+1.0e-4;
t1278 = t1248+1.0e-4;
t1279 = 1.0./sqrt(t1278);
t1280 = (p36.*t1279)./2.0;
t1281 = t1280-1.0./2.0;
t1282 = t1242+1.0e-4;
t1283 = 1.0./sqrt(t1282);
t1284 = (p33.*t1283)./2.0;
t1285 = t1284-1.0./2.0;
t1286 = pDOT32.*t1285.*5.0e1;
t1287 = -pDOT32+t609+t638+t842+t1066+t1286;
t1288 = L6.^2;
t1289 = (pDOT31.*t1253)./1.0e4;
t1290 = (pDOT32.*t1256)./1.0e4;
t1291 = (pDOT33.*t1259)./1.0e4;
t1292 = pDOT34.*t1262;
t1293 = pDOT35.*t1265;
t1294 = pDOT36.*t1268;
t1295 = t1243.^2;
t1296 = t1295./1.0e4;
t1297 = t1269.^2;
t1298 = t1297./1.0e4;
t1299 = t1270.^2;
t1300 = t1299./1.0e4;
t1323 = t1288./1.0e4;
t1324 = (pDOT34.*t1253)./1.0e4;
t1325 = (pDOT35.*t1256)./1.0e4;
t1326 = (pDOT36.*t1259)./1.0e4;
t1327 = pDOT31.*t1262;
t1328 = pDOT32.*t1265;
t1329 = pDOT33.*t1268;
t1301 = t1289+t1290+t1291+t1292+t1293+t1294+t1296+t1298+t1300-t1323-t1324-t1325-t1326-t1327-t1328-t1329;
t1302 = 1.0./sqrt(t1272);
t1303 = p36.*t1302.*2.5e4;
t1304 = t1303-1.0./2.0;
t1305 = pDOT36.*t1304.*3.0e2;
t1306 = sqrt(t1272);
t1307 = t1306./2.0;
t1330 = p36.*2.5e4;
t1308 = -pDOT36+t304+t314+t874+t1098+t1305+t1307-t1330-1.4715;
t1309 = 1.0./sqrt(t1277);
t1310 = p33.*t1309.*2.5e4;
t1311 = t1310-1.0./2.0;
t1312 = pDOT33.*t1311.*3.0e2;
t1313 = sqrt(t1277);
t1314 = t1313./2.0;
t1322 = p33.*2.5e4;
t1315 = -pDOT33+t648+t658+t873+t1097+t1312+t1314-t1322-1.4715;
t1316 = pDOT34.*t1281.*5.0e1;
t1317 = -pDOT34+t284+t294+t860+t1084+t1316;
t1318 = pDOT35.*t1281.*5.0e1;
t1319 = -pDOT35+t265+t289+t843+t1067+t1318;
t1320 = pDOT31.*t1285.*5.0e1;
t1321 = -pDOT31+t628+t633+t859+t1083+t1320;
t1331 = (t1243.*t1250.*t1253)./4.0;
t1332 = t1331-1.0;
t1333 = (t1243.*t1250.*t1301)./4.0;
t1334 = t1243.*t1250.*t1259.*t1308.*(5.0./3.0);
t1335 = t1243.*t1250.*t1256.*t1319.*(5.0./3.0);
t1336 = (t1250.*t1256.*t1269)./4.0;
t1337 = t1336-1.0;
t1338 = (t1250.*t1269.*t1301)./4.0;
t1339 = t1250.*t1259.*t1269.*t1308.*(5.0./3.0);
t1340 = t1250.*t1253.*t1269.*t1317.*(5.0./3.0);
t1341 = (t1250.*t1259.*t1270)./4.0;
t1342 = t1341-1.0;
t1343 = (t1250.*t1270.*t1301)./4.0;
t1344 = t1250.*t1253.*t1270.*t1317.*(5.0./3.0);
t1345 = t1250.*t1256.*t1270.*t1319.*(5.0./3.0);
out1 = [t333+t334+t335+t295.*t332.*(2.0e1./3.0)-t51.*t58.*t61.*t285.*(5.0./3.0)-t51.*t58.*t64.*t290.*(5.0./3.0)-t51.*t58.*t67.*t305.*(5.0./3.0);t338+t339+t340+t266.*t337.*(2.0e1./3.0)-t58.*t61.*t77.*t285.*(5.0./3.0)-t58.*t64.*t77.*t290.*(5.0./3.0)-t58.*t67.*t77.*t305.*(5.0./3.0);t343+t344+t345+t315.*t342.*(2.0e1./3.0)-t58.*t61.*t78.*t285.*(5.0./3.0)-t58.*t64.*t78.*t290.*(5.0./3.0)-t58.*t67.*t78.*t305.*(5.0./3.0);-t333-t334-t335+t285.*t332.*(2.0e1./3.0)+t51.*t58.*t64.*t290.*(5.0./3.0)-t51.*t58.*t61.*t295.*(5.0./3.0)+t51.*t58.*t67.*t305.*(5.0./3.0);-t338-t339-t340+t290.*t337.*(2.0e1./3.0)-t58.*t64.*t77.*t266.*(5.0./3.0)+t58.*t61.*t77.*t285.*(5.0./3.0)+t58.*t67.*t77.*t305.*(5.0./3.0);-t343-t344-t345+t305.*t342.*(2.0e1./3.0)+t58.*t61.*t78.*t285.*(5.0./3.0)+t58.*t64.*t78.*t290.*(5.0./3.0)-t58.*t67.*t78.*t315.*(5.0./3.0);t677+t678+t679+t629.*t676.*(2.0e1./3.0)-t395.*t402.*t405.*t634.*(5.0./3.0)-t395.*t402.*t408.*t639.*(5.0./3.0)-t395.*t402.*t411.*t659.*(5.0./3.0);t682+t683+t684+t610.*t681.*(2.0e1./3.0)-t402.*t405.*t421.*t634.*(5.0./3.0)-t402.*t408.*t421.*t639.*(5.0./3.0)-t402.*t411.*t421.*t659.*(5.0./3.0);t687+t688+t689+t649.*t686.*(2.0e1./3.0)-t402.*t405.*t422.*t634.*(5.0./3.0)-t402.*t408.*t422.*t639.*(5.0./3.0)-t402.*t411.*t422.*t659.*(5.0./3.0);-t677-t678-t679+t634.*t676.*(2.0e1./3.0)-t395.*t402.*t405.*t629.*(5.0./3.0)+t395.*t402.*t408.*t639.*(5.0./3.0)+t395.*t402.*t411.*t659.*(5.0./3.0);-t682-t683-t684+t639.*t681.*(2.0e1./3.0)-t402.*t408.*t421.*t610.*(5.0./3.0)+t402.*t405.*t421.*t634.*(5.0./3.0)+t402.*t411.*t421.*t659.*(5.0./3.0);-t687-t688-t689+t659.*t686.*(2.0e1./3.0)+t402.*t405.*t422.*t634.*(5.0./3.0)+t402.*t408.*t422.*t639.*(5.0./3.0)-t402.*t411.*t422.*t649.*(5.0./3.0);t901+t902+t903+t861.*t900.*(2.0e1./3.0)-t715.*t722.*t725.*t864.*(5.0./3.0)-t715.*t722.*t728.*t867.*(5.0./3.0)-t715.*t722.*t731.*t883.*(5.0./3.0);t906+t907+t908+t844.*t905.*(2.0e1./3.0)-t722.*t725.*t741.*t864.*(5.0./3.0)-t722.*t728.*t741.*t867.*(5.0./3.0)-t722.*t731.*t741.*t883.*(5.0./3.0);t911+t912+t913+t875.*t910.*(2.0e1./3.0)-t722.*t725.*t742.*t864.*(5.0./3.0)-t722.*t728.*t742.*t867.*(5.0./3.0)-t722.*t731.*t742.*t883.*(5.0./3.0);-t901-t902-t903+t864.*t900.*(2.0e1./3.0)-t715.*t722.*t725.*t861.*(5.0./3.0)+t715.*t722.*t728.*t867.*(5.0./3.0)+t715.*t722.*t731.*t883.*(5.0./3.0);-t906-t907-t908+t867.*t905.*(2.0e1./3.0)-t722.*t728.*t741.*t844.*(5.0./3.0)+t722.*t725.*t741.*t864.*(5.0./3.0)+t722.*t731.*t741.*t883.*(5.0./3.0);-t911-t912-t913+t883.*t910.*(2.0e1./3.0)+t722.*t725.*t742.*t864.*(5.0./3.0)+t722.*t728.*t742.*t867.*(5.0./3.0)-t722.*t731.*t742.*t875.*(5.0./3.0);t1125+t1126+t1127+t1085.*t1124.*(2.0e1./3.0)-t939.*t946.*t949.*t1088.*(5.0./3.0)-t939.*t946.*t952.*t1091.*(5.0./3.0)-t939.*t946.*t955.*t1107.*(5.0./3.0);t1130+t1131+t1132+t1068.*t1129.*(2.0e1./3.0)-t946.*t949.*t965.*t1088.*(5.0./3.0)-t946.*t952.*t965.*t1091.*(5.0./3.0)-t946.*t955.*t965.*t1107.*(5.0./3.0);t1135+t1136+t1137+t1099.*t1134.*(2.0e1./3.0)-t946.*t949.*t966.*t1088.*(5.0./3.0)-t946.*t952.*t966.*t1091.*(5.0./3.0)-t946.*t955.*t966.*t1107.*(5.0./3.0);-t1125-t1126-t1127+t1088.*t1124.*(2.0e1./3.0)-t939.*t946.*t949.*t1085.*(5.0./3.0)+t939.*t946.*t952.*t1091.*(5.0./3.0)+t939.*t946.*t955.*t1107.*(5.0./3.0);-t1130-t1131-t1132+t1091.*t1129.*(2.0e1./3.0)-t946.*t952.*t965.*t1068.*(5.0./3.0)+t946.*t949.*t965.*t1088.*(5.0./3.0)+t946.*t955.*t965.*t1107.*(5.0./3.0);-t1135-t1136-t1137+t1107.*t1134.*(2.0e1./3.0)+t946.*t949.*t966.*t1088.*(5.0./3.0)+t946.*t952.*t966.*t1091.*(5.0./3.0)-t946.*t955.*t966.*t1099.*(5.0./3.0);t1229+t1230+t1231-t1217.*t1228.*(2.0e1./3.0)-t1139.*t1146.*t1152.*t1183.*(5.0./3.0)+t1139.*t1146.*t1149.*t1213.*(5.0./3.0)-t1139.*t1146.*t1155.*t1211.*(5.0./3.0);t1234+t1235+t1236-t1183.*t1233.*(2.0e1./3.0)-t1146.*t1149.*t1165.*t1217.*(5.0./3.0)-t1146.*t1155.*t1165.*t1211.*(5.0./3.0)+t1146.*t1152.*t1165.*t1215.*(5.0./3.0);t1239+t1240+t1241-t1211.*t1238.*(2.0e1./3.0)-t1146.*t1152.*t1166.*t1183.*(5.0./3.0)+t1146.*t1155.*t1166.*t1204.*(5.0./3.0)-t1146.*t1149.*t1166.*t1217.*(5.0./3.0);-t1229-t1230-t1231-t1213.*t1228.*(2.0e1./3.0)+t1139.*t1146.*t1152.*t1183.*(5.0./3.0)+t1139.*t1146.*t1149.*t1217.*(5.0./3.0)+t1139.*t1146.*t1155.*t1211.*(5.0./3.0);-t1234-t1235-t1236-t1215.*t1233.*(2.0e1./3.0)+t1146.*t1152.*t1165.*t1183.*(5.0./3.0)+t1146.*t1149.*t1165.*t1217.*(5.0./3.0)+t1146.*t1155.*t1165.*t1211.*(5.0./3.0);-t1239-t1240-t1241-t1204.*t1238.*(2.0e1./3.0)+t1146.*t1152.*t1166.*t1183.*(5.0./3.0)+t1146.*t1149.*t1166.*t1217.*(5.0./3.0)+t1146.*t1155.*t1166.*t1211.*(5.0./3.0);t1333+t1334+t1335-t1321.*t1332.*(2.0e1./3.0)-t1243.*t1250.*t1256.*t1287.*(5.0./3.0)+t1243.*t1250.*t1253.*t1317.*(5.0./3.0)-t1243.*t1250.*t1259.*t1315.*(5.0./3.0);t1338+t1339+t1340-t1287.*t1337.*(2.0e1./3.0)-t1250.*t1253.*t1269.*t1321.*(5.0./3.0)-t1250.*t1259.*t1269.*t1315.*(5.0./3.0)+t1250.*t1256.*t1269.*t1319.*(5.0./3.0);t1343+t1344+t1345-t1315.*t1342.*(2.0e1./3.0)-t1250.*t1256.*t1270.*t1287.*(5.0./3.0)+t1250.*t1259.*t1270.*t1308.*(5.0./3.0)-t1250.*t1253.*t1270.*t1321.*(5.0./3.0);-t1333-t1334-t1335-t1317.*t1332.*(2.0e1./3.0)+t1243.*t1250.*t1256.*t1287.*(5.0./3.0)+t1243.*t1250.*t1253.*t1321.*(5.0./3.0)+t1243.*t1250.*t1259.*t1315.*(5.0./3.0);-t1338-t1339-t1340-t1319.*t1337.*(2.0e1./3.0)+t1250.*t1256.*t1269.*t1287.*(5.0./3.0)+t1250.*t1253.*t1269.*t1321.*(5.0./3.0)+t1250.*t1259.*t1269.*t1315.*(5.0./3.0);-t1343-t1344-t1345-t1308.*t1342.*(2.0e1./3.0)+t1250.*t1256.*t1270.*t1287.*(5.0./3.0)+t1250.*t1253.*t1270.*t1321.*(5.0./3.0)+t1250.*t1259.*t1270.*t1315.*(5.0./3.0)];
