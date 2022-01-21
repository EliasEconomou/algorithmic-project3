# -*- coding: utf-8 -*-
"""detect.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/17A0k-dW4B-ijsbDax9ZZ9dYC-8LgeIk0
"""



import sys
import math
import matplotlib.pyplot as plt
import keras
import pandas as pd
import numpy as np
import seaborn as sns
from keras.models import Sequential
from keras.layers import Dense
from keras.layers import LSTM
from keras.layers import Dropout
from keras.layers import *
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_absolute_error
from sklearn.model_selection import train_test_split
from keras.callbacks import EarlyStopping
pd.options.mode.chained_assignment = None

window = 30                 # default number of time steps when training
num_epochs = 6                # default number of epochs used when training the model                     
predict_series = 5         # default number of series to predict
dataset_path = "nasdaq2007_17.csv"
mae = 0.3 #default
for i in range(len(sys.argv)):
  if(sys.argv[i] == "-n"):
    predict_series = int(sys.argv[i+1])
  elif(sys.argv[i] == "-w"):
    window = int(sys.argv[i+1])
  elif(sys.argv[i] == "-e"):
    epochs = int(sys.argv[i+1])
  elif(sys.argv[i] == "-d"):
    dataset_path = sys.argv[i+1]
  elif(sys.argv[i] == "-mae"):
    dataset_path = sys.argv[i+1]

df=pd.read_csv(dataset_path, header=None, delimiter='\t') #create data frame from our csv file

df = df.transpose() # transpose rows to columns

# Renaming header as the time series' ids
df.columns = df.iloc[0]
df = df.reindex(df.index.drop(0)).reset_index(drop=True)
df.columns.name = None

num_series_selected = 100 # or df.shape[1] to select all series from dataset
print("Number of selected series to plot anomalies: ",num_series_selected)
series_to_train = []
for i in range(num_series_selected):
  series_to_train.append(i)
# print(series_to_train)

# Drop any columns (series) that have not been selected
df.drop(df.iloc[:, num_series_selected:], inplace = True, axis = 1)

# df

split = 0.8 # "split" percent of a series for training - the rest for testing
train_size = int(len(df) * split)
test_size = len(df) - train_size
train_series_list = []
test_series_list = []

# split every series into train/test part and apply standard scaler
# then concatenate them to one large series
for series_number in series_to_train:
  train_temp, test_temp = df.iloc[0:train_size,series_number:series_number+1], df.iloc[train_size:len(df),series_number:series_number+1]

  scaler = StandardScaler()
  scaler = scaler.fit(train_temp[[df.columns[series_number]]])
  train_temp[df.columns[series_number]] = scaler.transform(train_temp[[df.columns[series_number]]])
  test_temp[df.columns[series_number]] = scaler.transform(test_temp[[df.columns[series_number]]])

  train_series_list = train_series_list + train_temp[df.columns[series_number]].tolist()
  test_series_list = test_series_list + test_temp[df.columns[series_number]].tolist()

# create final dataframes train and test from lists
train = pd.DataFrame(train_series_list, columns =['value'])
test = pd.DataFrame(test_series_list, columns =['value'])

# Name 'time' the index column with incremented days
train.index.name = 'time'
test.index.name = 'time'

# function to reshape input
def create_dataset(X, y, time_steps=1):
    Xs, ys = [], []
    for i in range(len(X) - time_steps):
        v = X.iloc[i:(i + time_steps)].values
        Xs.append(v)
        ys.append(y.iloc[i + time_steps])
    return np.array(Xs), np.array(ys)

X_train, y_train = create_dataset(train[['value']],train.value,window)
X_test, y_test = create_dataset(test[['value']],test.value,window)

X_train.shape, y_train.shape, X_test.shape, y_test.shape

model = keras.Sequential()
model.add(keras.layers.LSTM(
    units=64,
    input_shape=(X_train.shape[1], X_train.shape[2])
))
model.add(keras.layers.Dropout(rate=0.2))
model.add(keras.layers.RepeatVector(n=X_train.shape[1]))
model.add(keras.layers.LSTM(units=64, return_sequences=True))
model.add(keras.layers.Dropout(rate=0.2))
model.add(keras.layers.TimeDistributed(keras.layers.Dense(units=X_train.shape[2])))
model.compile(loss='mae', optimizer='adam')

history = model.fit(
    X_train, y_train,
    epochs=num_epochs,
    batch_size=32,
    validation_split=0.1,
    shuffle=False
)

# model.save('detect_model.h5')

plt.plot(history.history['loss'],label='train')
plt.plot(history.history['val_loss'],label='test')
plt.legend()

X_train_pred = model.predict(X_train)
train_mae_loss = np.mean(np.abs(X_train_pred - X_train), axis=1)
# print(max(train_mae_loss))

sns.displot(train_mae_loss,bins=50,kde=True);

THRESHOLD = mae



print("Printing all n series that have anomalies:")
for series_number in range(20):
  test = df.iloc[train_size:len(df),series_number:series_number+1]
  test[df.columns[series_number]] = scaler.transform(test[[df.columns[series_number]]])
  
  test = test.rename(columns={test.columns[0]: 'value'})
  test.index.name = 'time'
  x_test, y_test = create_dataset(test[['value']],test.value,window)
  x_test_pred = model.predict(x_test)
  test_mae_loss = np.mean(np.abs(x_test_pred - x_test), axis=1)


  test_score_df = pd.DataFrame(index=test[window:].index)
  test_score_df['loss'] = test_mae_loss
  test_score_df['threshold'] = THRESHOLD
  test_score_df['anomaly'] = test_score_df.loss > test_score_df.threshold
  test_score_df['value'] = test[window:].value


  anomalies = test_score_df[test_score_df.anomaly == True]
  if anomalies.empty:
    continue


  value_df = test[window:]

  x = np.array(value_df)
  x.shape
  value_df = value_df.drop('value', 1)

  x = scaler.inverse_transform(x)

  value_df['value'] = x


  anom_df = anomalies

  anom_df = anom_df.drop('loss', 1)
  anom_df = anom_df.drop('threshold', 1)
  anom_df = anom_df.drop('anomaly', 1)

  a = np.array(anom_df)

  a = scaler.inverse_transform(a)
  anom_df['value'] = a


  plt.plot(
    test[window:].index, 
    value_df.value, 
    label='value'
  );
  # print(test[window:])
  sns.scatterplot(
    anomalies.index,
    anom_df.value,
    color=sns.color_palette()[3],
    s=52,
    label='anomaly'
  )
  plt.xticks(rotation=25)
  plt.legend()
  plt.show()

