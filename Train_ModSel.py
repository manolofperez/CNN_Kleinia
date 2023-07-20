import sys, os
import numpy as np
import time
import random
from random import shuffle, choice
import keras
import tensorflow as tf
from tensorflow.keras.models import Model
from tensorflow.keras.regularizers import l2
from tensorflow.keras.constraints import max_norm
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras import backend as K
from tensorflow.keras.layers import *
from tensorflow.keras.optimizers import *
from tensorflow.keras.models import load_model
from keras.utils import np_utils
import sklearn.metrics as metrics
from sklearn.metrics import log_loss
from sklearn.metrics import confusion_matrix

batch_size = 500

epochs = 500

num_classes = 4

#Define the CNN architecture
def create_cnn(xtest):
	inputShape = (xtest.shape[1], xtest.shape[2])
	inputs = Input(shape=inputShape)
	x = inputs
	x = Conv1D(125, kernel_size=3, activation='relu', use_bias=False, input_shape=(xtest.shape[1], xtest.shape[2]))(x)
	x = BatchNormalization()(x)
	x = Conv1D(250, kernel_size=3, use_bias=False, activation='relu')(x)
	x = BatchNormalization()(x)
	x = Conv1D(250, kernel_size=3, use_bias=False, activation='relu')(x)
	x = BatchNormalization()(x)
	x = MaxPooling1D(pool_size=3)(x)
	x = Flatten()(x)
	x = Dense(125, activation='relu')(x)
	x = Dropout(0.5)(x)
	x = Dense(125, activation='relu')(x)
	x = Dropout(0.5)(x)
	# The final fully-connected layer head will have a softmax dense layer
	logits = Dense(num_classes, name='logits')(x)
	out = Activation('softmax')(logits)

	# Construct the CNN
	model = Model(inputs, out)
	# Return the CNN
	return model

SS = np.load("../testSet/Mod_SS.npy",mmap_mode='r')
SSHclim = np.load("../testSet/Mod_SSHclim.npy",mmap_mode='r')
CIH1 = np.load("../testSet/Mod_CIH1.npy",mmap_mode='r')
CIH2 = np.load("../testSet/Mod_CIH2.npy",mmap_mode='r')

X=np.concatenate((SS['Mod_SS'],SSHclim['Mod_SSHclim'],CIH1['Mod_SSHclim'],CIH2['Mod_SSHclim']),axis=0)

#transform major alleles in -1 and minor 1
for arr,array in enumerate(X):
  for idx,row in enumerate(array):
    if np.count_nonzero(row) > len(row)/2:
      X[arr][idx][X[arr][idx] == 1] = -1
      X[arr][idx][X[arr][idx] == 0] = 1
    else:
      X[arr][idx][X[arr][idx] == 0] = -1


y=[0 for i in range(len(SS['Mod_SS']))]
y.extend([1 for i in range(len(SSHclim['Mod_SSHclim']))])
y.extend([2 for i in range(len(CIH1['Mod_SSHclim']))])
y.extend([3 for i in range(len(CIH2['Mod_SSHclim']))])

y = np.array(y)

del (SS,SSHclim,CIH1,CIH2)
print (len(X), len(y))
shf = list(range(len(X)))
shuffle(shf)

y_shf = y[shf]
X = X[shf]

#Add missing data as 0s, according to a specifies missing data percentage
#215,160 SNPs and 2,562 missing genotypes = 4.4%
missD_perc = 1.19
missD = int(X.shape[1]*X.shape[2]*(missD_perc/100))
for i in range(X.shape[0]):
  for m in range(missD):
    j = random.randint(0, X.shape[1] - 1)
    k = random.randint(0, X.shape[2] - 1)
    X[i][j][k] = 0
del(missD)

X=np.array(X)

xtrain, xtest = X[int(len(y)*.25):], X[:int(len(y)*.25)]
ytrain, ytest = y_shf[int(len(y_shf)*.25):], y_shf[:int(len(y_shf)*.25)]

ytest = np_utils.to_categorical(ytest, num_classes)
ytrain = np_utils.to_categorical(ytrain, num_classes)

# Create the CNN network
model = create_cnn(xtest)

opt = SGD(learning_rate=0.001)

model.compile(loss=keras.losses.categorical_crossentropy,
	              optimizer=opt,
	              metrics=['accuracy'])

model.summary()

earlyStopping = EarlyStopping(monitor='val_accuracy', patience=150, verbose=0, mode='max', restore_best_weights=True)
start = time.time()
model.fit(xtrain, ytrain, batch_size=batch_size,
          epochs=epochs,
          verbose=1,
          validation_data=(xtest, ytest),callbacks=[earlyStopping])
print (f'Time: {time.time() - start}')

model.save(filepath='Trained_Kleinia.acc.mod')

pred = model.predict(xpred)
pred_cat = [i.argmax() for i in pred]
print (confusion_matrix(y, pred_cat))
print (confusion_matrix(y, pred_cat) / float(len(y)))

np.savetxt("testSet_Predictions.txt", pred)

infile=np.loadtxt("../input_Kleinia.txt")
inp=np.array(infile)
num_samples=100
res = []
for i in range(0,num_samples):
	idx = np.random.choice(inp.shape[0], 1000, replace=False)
	n = inp[idx,:]
	res.append(np.array(n))

Emp = np.array(res)
Emp_pred = model.predict(Emp)

np.savetxt("Emp_Predictions.txt", Emp_pred)