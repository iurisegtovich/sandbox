import os,sys
import numpy as np
import argparse


def knnC(Dtrain,Dtest,k):
    prediction = -1*np.ones(len(Dtest))
    success = 0
    for i in range(len(Dtest)):
        distance = np.sqrt(((Dtrain[:,:-1]-Dtest[i,:-1])**2).sum(1))
        I = np.argsort(distance)[0:k]
        if Dtrain[I,-1].mean() >= 0.5:
            prediction[i] =1
        else:
            prediction[i] =0

        if prediction[i] == Dtest[i,-1]:
            success += 1.0

    return prediction,success/len(Dtest)

def knnR(Dtrain,Dtest,k):
    prediction = -1*np.ones(len(Dtest))
    for i in range(len(Dtest)):
        distance = np.sqrt(((Dtrain[:,:-1]-Dtest[i,:-1])**2).sum(1))
        I = np.argsort(distance)[0:k]
        prediction[i] = Dtrain[I,-1].mean()
    MSE = ((prediction - Dtest[:,-1])**2).mean()

    return prediction,MSE

def lr_learn(X,Y):
    # add dummuy variable x0=1
    X1 = np.ones((X.shape[0],X.shape[1]+1))
    X1[:,1:] = X[:,:]
    X1 = np.matrix(X1)
    Y = np.matrix(Y)
    Y = Y.reshape(-1,1)
    model = (X1.T*X1).I*X1.T*Y
    return model

def lr_predict(X,Y,model):
    X1 = np.ones((X.shape[0],X.shape[1]+1))
    X1[:,1:] = X[:,:]
    X1 = np.matrix(X1)
    prediction = X1*model
    prediction = np.array(prediction).reshape(1,-1)
    MSE = ((prediction - Y)**2).mean()
    return prediction,MSE




parser = argparse.ArgumentParser()
parser.add_argument('-tr','--Training',help='Path to the training data',required=True)
parser.add_argument('-tt','--Testing',help='Path to the testing data',required=True)
algorithm_parser = parser.add_subparsers()

parser_kNNC = algorithm_parser.add_parser('kNNC',help='Use k-nearest neighbor classification algorithm')
parser_kNNC.add_argument('-k',type=int,help='k-nearest neighbours')
parser_kNNC.set_defaults(func='knnC')

parser_kNNR = algorithm_parser.add_parser('kNNR',help='Use k-nearest neighbor regression algorithm')
parser_kNNR.add_argument('-k',type=int,help='k-nearest neighbours')
parser_kNNR.set_defaults(func='knnR')

parser_LR = algorithm_parser.add_parser('LR',help='Use linear regression algorithm')
parser_LR.add_argument('-n',type=int,help='number of varialbes')
parser_LR.set_defaults(func='LR')

args=parser.parse_args()

Dtrain = np.loadtxt(args.Training)
Dtest = np.loadtxt(args.Testing)

if args.func == 'knnC':
    prediction, accuracy = knnC(Dtrain,Dtest,args.k)
    print prediction,accuracy

if args.func == 'knnR':
    prediction, MSE = knnR(Dtrain,Dtest,args.k)
    print prediction,MSE

if args.func == 'LR':
    Xtrain = Dtrain[:,:-1]
    Ytrain = Dtrain[:,-1]
    Xtest = Dtest[:,:-1]
    Ytest = Dtest[:,-1]
    model = lr_learn(Xtrain,Ytrain)
    prediction,MSE = lr_predict(Xtest,Ytest,model)
    print prediction,MSE


