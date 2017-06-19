function [ flag ] = predicate_OWN(region1, region2)
    sd1 = std2(region1(region1 ~= 0));
    m1 = mean2(region1(region1 ~= 0));
    sd2 = std2(region2(region2 ~= 0));
    m2 = mean2(region2(region2 ~= 0));
    if (m1 ~= m1) m1 = 0; end % dealing with NaN s from all-zeros regions
    if (m2 ~= m2) m2 = 0; end

    flag = (abs(m2-m1) < 40); % parametr

end
