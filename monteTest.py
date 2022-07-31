import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv('monte data 10000.csv',sep = ' ')
sns.distplot(df['a'])
plt.show()