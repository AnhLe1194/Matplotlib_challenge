#!/usr/bin/env python
# coding: utf-8

# ## Observations and Insights
# #### Add your analysis here
# ---

# In[2]:


# Dependencies and Setup
import matplotlib.pyplot as plt
from scipy.stats import linregress
import pandas as pd
import scipy.stats as st
import numpy as np

get_ipython().run_line_magic('matplotlib', 'notebook')
# Study data files
mouse_metadata_path = "data/Mouse_metadata.csv"
study_results_path = "data/Study_results.csv"

# Read the mouse data and the study results
mouse_metadata = pd.read_csv(mouse_metadata_path)
study_results = pd.read_csv(study_results_path)

# Combine the data into a single dataset
merge_df = pd.merge(mouse_metadata, study_results, on= "Mouse ID", how= "outer")
# Display the data table for preview
mouse_metadata


# In[3]:


study_results


# In[4]:


merge_df


# In[5]:


# Check the number of mice.
mice = mouse_metadata['Mouse ID'].value_counts()
number_of_mice = len(mice)
print(number_of_mice)


# In[6]:


# Getting the duplicate mice by ID number that shows up for Mouse ID and Timepoint.
dup_mouse_id = merge_df.loc[merge_df.duplicated(subset=['Mouse ID', 'Timepoint',]),'Mouse ID'].unique()
dup_mouse_id


# In[7]:


dup = merge_df["Mouse ID"] == 'g989'


# In[8]:


# Optional: Get all the data for the duplicate mouse ID. 


# In[17]:


# Create a clean DataFrame by dropping the duplicate mouse by its ID.
clean_data=merge_df[merge_df["Mouse ID"].isin(dup_mouse_id) == False]
clean_data


# In[18]:


# Check the number of mice in the clean DataFrame.
mice = clean_data['Mouse ID'].unique() 
number_of_mice = len(mice)
print(number_of_mice)


# In[19]:


mice
mice_id = pd.DataFrame


# ## Summary Statistics

# In[20]:


# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen
grouped_drug = merge_df.groupby('Drug Regimen')


# Use this straighforward method, create multiple series and put them all in a dataframe at the end.
grouped_drug = grouped_drug.agg(['mean','median','var','std','sem'])["Tumor Volume (mm3)"]
grouped_drug


# In[21]:


# Generate a summary statistics table of mean, median, variance, standard deviation, and SEM of the tumor volume for each regimen

grouped_drug = merge_df.groupby(['Drug Regimen'])

mean = grouped_drug["Tumor Volume (mm3)"].mean()
mean                               

median = grouped_drug["Tumor Volume (mm3)"].median()
median

var = grouped_drug["Tumor Volume (mm3)"].var()

std = grouped_drug["Tumor Volume (mm3)"].std()

sem = grouped_drug["Tumor Volume (mm3)"].sem()

# Use method to produce everything with a single groupby function
sum = pd.DataFrame({"mean": mean, "median": median, "var": var, "std": std, "sem": sem})
sum


# In[33]:


plt.figure(1)
clean_data['Drug Regimen'].value_counts().plot(kind='bar')


# In[23]:


# Generate a bar plot showing the total number of mice for each treatment throughout the course of the study using pyplot.
treatment = clean_data["Drug Regimen"].value_counts()

# Set x axis and tick locations
x_axis = np.arange(len(treatment))
tick_locations = [value+0.4 for value in x_axis]

# Create a list indicating where to write x labels and set figure size to adjust for space
plt.figure(2)
plt.bar(x_axis, treatment, color='b', alpha=0.25, align="edge")
plt.xticks(tick_locations, treatment.index, rotation="45")

# Set x and y limits
plt.xlim(-0.5, len(x_axis))
plt.ylim(0, max(treatment)+10)

plt.title("Total number of mice")
plt.xlabel("Drug Regimen")
plt.ylabel("Count")


# In[24]:


# Generate a pie plot showing the distribution of female versus male mice using pandas

list = {"Sex": ["Male", "Female"], "Percentage": [0.45, 0.55]}
df = pd.DataFrame(list)

Sex = clean_data.groupby(['Sex'])
Sex['Sex'].value_counts()


# In[25]:


only_male = clean_data.loc[clean_data["Sex"]== "Male", :]
male = only_male["Mouse ID"].unique()
total_male = len(male)
per_male = total_male*100/number_of_mice

only_female = clean_data.loc[clean_data["Sex"]== "Female", :]
female = only_female["Mouse ID"].unique()
total_female = len(female)
per_female = total_female*100/number_of_mice

gender_df = pd.DataFrame([[total_female, per_female], [total_male, per_male]], index=["Female","Male"], columns=["Total Count","Percentage of Mouse"])
gender_df


# In[26]:


colors = ['blue', 'red']
explode = (0.1, 0)
plot = gender_df.plot.pie(y='Total Count',figsize=(5,5), colors = colors, startangle=140, explode = explode, shadow = True, autopct="%1.1f%%")


# In[34]:


# Generate a pie plot showing the distribution of female versus male mice using pyplot
# Labels for the sections of our pie chart
labels = ["Female", "Male"]

# The values of each section of the pie chart
sizes = [49.6, 50.4]

# The colors of each section of the pie chart
colors = ["red", "lightskyblue"]

# Tells matplotlib to seperate the "Humans" section from the others
explode = (0.1, 0)

plt.figure(4)
plt.pie(sizes, explode=explode, labels=labels, colors=colors,
        autopct="%1.1f%%", shadow=True, startangle=140)

# Create axes which are equal so we have a perfect circle
plt.axis("equal")
plt.show()


# 
# ## Quartiles, Outliers and Boxplots

# In[35]:


# Calculate the final tumor volume of each mouse across each of the treatment regimens:
# Start by getting the last (greatest) timepoint for each mouse
# Merge this group df with the original dataframe to get the tumor volume at the last timepoint


# In[36]:


mouse_id = clean_data.groupby(["Mouse ID", "Drug Regimen"])
mouse_id.head()
max_timepoint = mouse_id["Timepoint"].max()
max_timepoint

df_1 = pd.DataFrame({"Timepoint": max_timepoint})
df_1

df_1 = df_1.reset_index()
df_1

last_tumor = pd.merge(df_1, clean_data, on= ("Mouse ID", "Timepoint", "Drug Regimen"), how= "left")
last_tumor


# In[37]:


only_Capomulin = merge_df.loc[merge_df['Drug Regimen'] == "Capomulin", :]
only_Capomulin
Capomulin = only_Capomulin['Tumor Volume (mm3)']
only_Capomulin


# In[55]:


# Put 4 treatment names into a list for use with a for loop (and later for plot labels)
treatment_list = ["Capomulin", "Ramicane", "Infubinol", "Ceftamin"]

# Create a empty list to fill with tumor vol data (for plotting) (hint: each element of the list will be series)
tumor_vol_list = []

# For each treatment in the list, calculate the IQR and quantitatively 
# determine if there are any potential outliers. 
# Locate the rows which contain mice on each drug and get the tumor volumes
# Determine outliers using upper and lower bounds

for drug in treatment_list:
    only_drug = merge_df.loc[merge_df['Drug Regimen'] == drug, :]
    drug_name = only_drug["Tumor Volume (mm3)"]
    quartiles = drug_name.quantile([.25,.5,.75])
    lowerq = quartiles[0.25]
    upperq = quartiles[0.75]
    iqr = upperq-lowerq
    
    print(f"The lower quartile of tumor sizes is: {lowerq}")
    print(f"The upper quartile of tumor sizes is: {upperq}")
    print(f"The interquartile range of tumor sizes is: {iqr}")
    print(f"The the median of tumor sizes is: {quartiles[0.5]} ")
    
    lower_bound = lowerq - (1.5*iqr)
    upper_bound = upperq + (1.5*iqr)
    
    print(f"{drug}: Values below {lower_bound} could be outliers.")
    print(f"{drug}: Values above {upper_bound} could be outliers.")
    
    tumor_df = pd.DataFrame({"Tumor Volume (mm3)": drug_name})
    tumor_df
    
    tumor = pd.DataFrame({"Tumor Volume (mm3)": drug_name})
    tumor
    outliers = tumor.loc[(tumor["Tumor Volume (mm3)"] > upper_bound) | (tumor["Tumor Volume (mm3)"] < lower_bound ), :] 
    outliers
    print(f" Outliers: {outliers['Tumor Volume (mm3)']}")
    print("\n")


# In[51]:


only_drug = merge_df.loc[merge_df['Drug Regimen'] == "Infubinol", :]
drug_name = only_drug["Tumor Volume (mm3)"]
tumor_df = pd.DataFrame({"Tumor Volume (mm3)": drug_name})
tumor_df
outliers = tumor_df.loc[(tumor_df["Tumor Volume (mm3)"] > 72.31757996875001 ) | (tumor_df["Tumor Volume (mm3)"] < 32.309217298749985 ), :] 
outliers


# In[23]:


# Generate a box plot of the final tumor volume of each mouse across four regimens of interest
last_tumor


# In[54]:


tumors1 = last_tumor.loc[last_tumor['Drug Regimen'] == "Capomulin", :]
tumors_1 = tumors1["Tumor Volume (mm3)"]
tumors2 = last_tumor.loc[last_tumor['Drug Regimen'] == "Ramicane", :]
tumors_2 = tumors2["Tumor Volume (mm3)"]
tumors3 = last_tumor.loc[last_tumor['Drug Regimen'] == "Infubinol", :]
tumors_3 = tumors3["Tumor Volume (mm3)"]
tumors4 = last_tumor.loc[last_tumor['Drug Regimen'] == "Ceftamin", :]
tumors_4 = tumors4["Tumor Volume (mm3)"]

data_to_plot = [tumors_1, tumors_2, tumors_3, tumors_4]

plt.figure(5)
fig1, ax1 = plt.subplots()
ax1.set_title('Tumors')
ax1.set_ylabel('Final Tumor Volume (mm3)')
ax1.set_xlabel('Drug Regimen')

ax1.boxplot(data_to_plot, labels=["Capomulin","Ramicane","Infubinol","Ceftamin",])

plt.savefig('boxplot')
plt.show()


# ## Line and Scatter Plots

# In[53]:


# Generate a line plot of time point versus tumor volume for a mouse treated with Capomulin
merge_df = pd.merge(mouse_metadata, study_results, on= "Mouse ID", how= "left")
Capomulin = merge_df.loc[merge_df['Drug Regimen'] == "Capomulin", :]
tumor_vol = Capomulin.loc[Capomulin["Mouse ID"] == "l509",:]

x_axis = tumor_vol["Timepoint"]
y_axis = tumor_vol["Tumor Volume (mm3)"]

plt.figure(6)
plt.title('Capomulin treatmeant of mouse l509')
plt.plot(x_axis, y_axis,linewidth=2, markersize=12)
plt.xlabel('Timepoint (Days)')
plt.ylabel('Tumor Volume (mm3)')

plt.savefig('linechart')
plt.show()


# In[52]:


# Generate a scatter plot of mouse weight versus average tumor volume for the Capomulin regimen
avg = Capomulin.groupby(['Mouse ID']).mean()

plt.figure(7)
plt.scatter(avg['Weight (g)'], avg['Tumor Volume (mm3)'])
plt.xlabel('Weight (g)')
plt.ylabel('Average Tumor Volume (mm3)')

plt.savefig('scatterplot')
plt.show()


# ## Correlation and Regression

# In[55]:


# Calculate the correlation coefficient and linear regression model 
# for mouse weight and average tumor volume for the Capomulin regimen
correlation =round(st.pearsonr(avg['Weight (g)'],avg['Tumor Volume (mm3)'])[0],2)
print(f"correlation: {correlation}")


# In[56]:


model=st.linregress(avg['Weight (g)'],avg['Tumor Volume (mm3)'])
model


# In[59]:


x_values = avg['Weight (g)']
y_values = avg['Tumor Volume (mm3)']
(slope, intercept, rvalue, pvalue, stderr) = linregress(x_values, y_values)
regress_values = x_values * slope + intercept
line_eq = "y = " + str(round(slope,2)) + "x + " + str(round(intercept,2))

plt.figure(8)
plt.scatter(x_values,y_values)
plt.plot(x_values,regress_values,"r-")
plt.annotate(line_eq,(6,10),fontsize=15,color="red")
plt.xlabel('Weight(g)')
plt.ylabel('Average Tumore Volume (mm3)')
plt.show()

