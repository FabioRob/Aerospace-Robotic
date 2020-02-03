function [] = plotPlanetPosition (assumption ,planet_id ,date , marker , markersize , col )
% This function plot position of planet at given date

[~, r, ~, ~]= planet_elements_and_sv(assumption,planet_id, date(1),date(2),date(3),date(4),date(5),date(6));    
 plot3(r(1),r(2),r(3),'Marker',marker,'MarkerSize' ,markersize , 'color',col)
 