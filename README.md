# Layers
%xlsreader for the data to be implemented into LineThruPoly. Now you need
%a safe place to store your database in a csv file or something.

%first, sort by layers
layers_to_consider = xlsread('Geographical Considerations in Transmission Line Location','A:A');
[~,polygon_names] = xlsread('Geographical Considerations in Transmission Line Location', 'B:B');
layers = xlsread('Geographical Considerations in Transmission Line Location','C:C');
polygon_cost = xlsread('Geographical Considerations in Transmission Line Location','D:D');
tolerance = xlsread('Geographical Considerations in Transmission Line Location','D:D');
%l1 = length(polygon_names);
%l2 = length(layers);
%l3 = length(polygon_cost);
%l4 = length(tolerance);
%ls = {l4,l2,l3,l1};
%comp_matrix = eye(length(ls));
%for i = 1:length(ls)
 %   if i == 4
  %      if ls{i} == ls{1} +1
   %         comp_matrix(i,1) =1;
    %        comp_matrix(1,i) = 1;
     %   end
   % else
  %     if ls{i} == ls{i+1}
  %         comp_matrix(i,i+1) = 1;
   %         comp_matrix(i+1,i) = 1;
   %     end
   % end
%end
%if ~isempty(find(~comp_matrix,1))
%    error('input data missing')
%end
split = find(polygon_names == "Obstacle Name");
if split == 1 %if split is one then there is no case study shape boundary is submitted
    within_polygon = [];
else
    within_polygon  = polygon_names(2:split(1)-2);
end
if length(layers)>=split  %if there are blocking polygons to consider and an appropriate splittage to be made
    layers_within = layers(1:split(1)-3);
    layers_outside = layers(split:end);
    polygon_cost_within = polygon_cost(1:split(1)-3);
    polygon_cost_outside = polygon_cost(split:end);
    tolerance_within = tolerance(1:split(1)-3);
    tolerance_outside = tolerance(split:end);
elseif split == 1 %if split == 1 then no case study shape boundary was submitted
    layers_outside = layers;
    polygon_cost_outside = polygon_cost;
    tolerance_outside = tolerance;
else 
    layers_within = layers;
    polygon_cost_within = polygon_cost;
    tolerance_within = tolerance;
end
if split == length(polygon_names) %if no blocking polys are submitted
    outside_polygon = []; %it's empty and we won't do anything widit
else
    outside_polygon = polygon_names(split+1:end);
end
i = 1;
while i <= length(within_polygon)
    if within_polygon(i) ~= ""
        if ismember(layers_within(i), layers_to_consider) %if we consider it, then we don't remove it
            i = i + 1;
        else
            within_polygon(i) = [];
            if ~isempty(outside_polygon)
                layers_within(i) = [];
                polygon_cost_within(i) = [];
                tolerance_within(i) = [];
            end
        end
    else
        within_polygon(i) = [];
        if ~isempty(outside_polygon)
            layers_within(i) = [];
            polygon_cost_within(i) = [];
            tolerance_within(i) = [];
        end
    end
    
end
i = 1;
while i <= length(outside_polygon)
    if outside_polygon(i) ~= ""
        if ismember(layers_outside(i),layers_to_consider)
            i = i + 1;
        else
            outside_polygon(i) = [];
            layers_outside(i) = [];
            polygon_cost_outside(i) = [];
            tolerance_outside(i) = [];
        end
    else
        outside_polygon(i) = [];
        layers_outside(i) = [];
        polygon_cost_outside(i) = [];
        tolerance_outside(i) = [];
    end
end
%now call the necessary csvs
%polygonx = nan(length(within_polygon)+length(outside_polygon),100000);
%polygony = nan(length(within_polygon)+length(outside_polygon),100000);
%j = 1;
i = 1;
while i <= length(within_polygon) %for the within polygons, we want to add up all the strings and turn it into one big thing :) 
    polygon = csvread(within_polygon(i) + ".csv");
    country_vectorx_append = [polygon(:,1); nan];
    country_vectory_append = [polygon(:,2); nan];
    if i ==1 %in the first iteration we initialize the first polygon
       country_vectorx = country_vectorx_append;
       country_vectory = country_vectory_append;
    else %we concatenate the strings
        
        country_vectorx = [country_vectorx; country_vectorx_append];
        country_vectory = [country_vectory; country_vectory_append];
        %polygonx(j,find(isnan(polygonx(j,:)),1):size(polygon,1)) = [polygonx(j,:) polygon(:,1)'];
        %polygony(j,find(isnan(polygony(j,:)),1):size(polygon,1)) = [polygony(j,:) polygon(:,2)'];
    end
    if i == length(within_polygon)
        begindex(1) = 1;
        
    end
    i = i+1;
end
polygonx = country_vectorx;
polygony = country_vectory;
max_x = max(polygonx);
min_x = min(polygonx);
max_y = max(polygony);
min_y = min(polygony);
new_shape_x = [max_x+1, max_x+1, min_x - 1, min_x-1, max_x+1]';
new_shape_y = [min_y-1, max_y+1, max_y + 1, min_y-1, min_y-1]';
polygonx = [flipud(polygonx); nan; new_shape_x; nan];
polygony = [flipud(polygony); nan; new_shape_y; nan];
endex(1) = length(polygonx);
fclose('all');
i = 1;
while i <= length(outside_polygon) %for the outside polygons, we store the information for the blocking polygon
    polygon = csvread(outside_polygon(i) + ".csv");
    begindex(i+1) = length(polygonx)+1;
    polygonx = [polygonx; polygon(:,1)];
    polygony =  [polygony; polygon(:,2)];
    endex(i+1) = length(polygonx);
    i = i+1;
end
%replace all the zeros from expanding the matrix with nans
    lookfor = find(isnan(polygonx));
    j = length(within_polygon) + 2; %ignore the nans between country boundaries
    old_begindex = begindex;
    old_endex = endex;
while j <= length(lookfor)
    if ismember(lookfor(j)-1,lookfor) %if it's a repeat (oops)
        lookfor(j) = [];
    else
        begindex(i+1) = lookfor(j)+2;
        endex(i+1) = lookfor(j)-1;
        j = j + 1;
        i = i + 1;
    end
end
begindex = sort(begindex);
endex = sort(endex);
j = length(tolerance_outside)+1;
old_tolerance_outside = tolerance_outside;
old_polygon_cost_outside = polygon_cost_outside;
for i = 2:length(old_begindex)
    lookforbegin = find(begindex>=old_begindex(i));
    lookforend = find(endex<=old_endex(i));
    lookfor = intersect(lookforbegin,lookforend);
    new_tolerance_outside = old_tolerance_outside(i-1)*ones(length(lookfor),1);
    new_polygon_cost_outside = old_polygon_cost_outside(i-1)*ones(length(lookfor),1);
    if i == 2
        tolerance_outside = new_tolerance_outside;
        polygon_cost_outside = new_polygon_cost_outside;
    else
        tolerance_outside = [tolerance_outside; new_tolerance_outside];
        polygon_cost_outside = [polygon_cost_outside; new_polygon_cost_outside];
    end
end
mesh_spacing = xlsread('Geographical Considerations in Transmission Line Location','G4:G4');
lake_jump = xlsread('Geographical Considerations in Transmission Line Location','H4:H4');
proximity = xlsread('Geographical Considerations in Transmission Line Location','I4:I4');
shortestpath_bool = xlsread('Geographical Considerations in Transmission Line Location','J4:J4');
shortestpath_lb = xlsread('Geographical Considerations in Transmission Line Location', 'K4:K4');
shortestpath_ub = xlsread('Geographical Considerations in Transmission Line Location.xlsx','K5:K5');
reference_latitude = xlsread('Geographical Considerations in Transmission Line Location', 'L4:L4');
node_from_index_case_study = xlsread('Geographical Considerations in Transmission Line Location','NodeInfo','D:D');
node_to_index_case_study = xlsread('Geographical Considerations in Transmission Line Location','NodeInfo','H:H');
node_from_location_case_study = xlsread('Geographical Considerations in Transmission Line Location','NodeInfo','E:F');
node_to_location_case_study = xlsread('Geographical Considerations in Transmission Line Location','NodeInfo','I:J');

node_from_location_case_study = fliplr(node_from_location_case_study); %flip longtitude and latitudes
node_to_location_case_study = fliplr(node_to_location_case_study);

%define the outer boundary as bad and rest as inner boundary as good
boundary = zeros(length(begindex),1);
boundary(1) = 1;

%consider reference latitude. If nothing is entered then use mean radius
if isempty(reference_latitude)
    reference_latitude = 6371;
end
% consider mesh spacing. if empty, default to 1km
if isempty(mesh_spacing)
    mesh_spacing = 1;
end
% consider lake jump. if empty, default to 0. 
if isempty(lake_jump)
    lake_jump = 0;
end
% proximity to boundary check. if empty, default to 10^6km
if isempty(proximity)
    proximity = 1e+6;
end
% shortest path calculation if empty default to 1.
if isempty(shortestpath_bool)
    shortestpath_bool = 0;
end
% shortest path lower bound. if empty default to 0
if isempty(shortestpath_lb)
    shortestpath_lb = 0;
end
%shortest path upper bound. if empty default to 1e+6
if isempty(shortestpath_ub)
    shortestpath_ub = 1e+6;
end

polygons = [polygonx, polygony];
polygon_cost_within = mean(polygon_cost_within);
tolerance_within = mean(tolerance_within);
polycost_case_study = [polygon_cost_within; polygon_cost_outside];
tolerance = [tolerance_within; tolerance_outside];
[minimized_distance, straight_line_distance, penalized_straight_line_distance] = ...
    LineThruPoly2(polygons, begindex, endex, node_from_index_case_study, ...
    node_to_index_case_study, node_from_location_case_study , node_to_location_case_study, ... 
    polycost_case_study, mesh_spacing, boundary, tolerance, lake_jump, proximity, ...
    shortestpath_bool, shortestpath_lb, shortestpath_ub,reference_latitude);

table('Minimized Distance', 'Straight Line Distance', 'Penalized Straight Line Distance')
