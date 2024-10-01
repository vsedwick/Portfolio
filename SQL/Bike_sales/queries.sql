--query1: What bike model is most preferred by females compared to males on the west coast

select Bike_Model,
sum(case when Customer_Gender = 'Female' then 1 else 0 end) as Female,
sum(case when Customer_Gender = 'Male' then 1 else 0 end) as Male
from bike_sales
where Store_Location in ("San Antonio", "Houston", "Pheonix", "Los Angeles")
group by Bike_Model
order by Bike_Model;


--query2: which city has the highest number of bike purchases?

select 'which city has the highest number of bike purchases';

select Store_Location, sum(Quantity) as Total_Purchases
from bike_sales
group by Store_Location
order by Total_Purchases DESC
limit 1;


--query 3: What is the average bike price for different age groups
select 'What is the average bike price for different age groups';
select case when Customer_age between 18 and 25 then 'Young Adult (18-25)'
when Customer_age between 26 and 39 then 'Adult (26-39)'
when Customer_age between 40 and 57 then 'Middle Aged (40-57)'
when Customer_age between 58 and 70 then 'Elderly (58-70)'
else 'Other'
end as Age_groups,
avg(Price) as Averaged_by_age
from bike_sales
group by Age_groups
order by Averaged_by_age desc;


--query 4

select 'What is the sex split for different age groups';
select case when Customer_age between 18 and 25 then 'Young Adult (18-25)'
when Customer_age between 26 and 39 then 'Adult (26-39)'
when Customer_age between 40 and 57 then 'Middle Aged (40-57)'
when Customer_age between 58 and 70 then 'Elderly (58-70)'
else 'Other'
end as Age_groups,
sum(case when Customer_Gender = 'Female' then 1 else 0 end) as Female,
sum(case when Customer_Gender = 'Male' then 1 else 0 end) as Male
from bike_sales
group by Age_groups
order by Female desc;

--query 5: What payment method is most commonly used by younger vs. older customers?

select '--query 5: What payment method is most commonly used by younger vs. older customers (raw numbers)';

with Age_Group as (
    select case when Customer_age between 18 and 25 then 'Young Adult (18-25)'
    when Customer_age between 26 and 39 then 'Adult (26-39)'
    when Customer_age between 40 and 57 then 'Middle Aged (40-57)'
    when Customer_age between 58 and 70 then 'Elderly (58-70)'
    else 'Other'
    end as Age_groups,
    Payment_Method,
    count(*) as Method_Use_no
    from bike_sales
    group by Age_groups, Payment_Method
)
Select Age_groups, Payment_Method
from (
    select Age_groups, Payment_Method, rank() over (partition by 
    Age_groups order by  Method_Use_no desc) as rank
    from Age_Group
) as Ranked_payment_method
where rank = 1;


--query 6: Which bike model generates the most revenue?
select Bike_Model, sum(Price*Quantity) as Total_Revenue from bike_sales
group by Bike_Model
order by Total_Revenue desc
limit 1;

--query 7: Do customers of a specific gender tend to purchase more expensive bikes?

select Customer_Gender, avg(Price)
from bike_sales
group by Customer_Gender;

--query 8: Which store location has the highest revenue?

select Store_Location, sum(Price*Quantity) as Total_Revenue from bike_sales
group by Store_Location
order by Total_Revenue desc
limit 1;



