# %%
import pandas as pd
import numpy
from lightgbm.sklearn import LGBMClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix, accuracy_score, roc_curve, auc
from sklearn.ensemble import RandomForestClassifier
import catboost
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.naive_bayes import BernoulliNB, GaussianNB, MultinomialNB


# %%
# model
def model_algorithm(file_data, file_group, test_size_tmp, model):
    data_t = file_data.T
    data_group = file_group
    
    # X = data_t.iloc[:, :-1].values
    # Y = data_t.iloc[:, -1].values

    # X = data_t[:, :-1]
    # Y = data_t[:, -1]
    
    # 创建一个映射字典，将'Normal'映射为1，'Tumor'映射为2
    # replacement_dict = {'normal': 0, 'tumor': 1}

    # 使用.replace()方法替换'Sample'列中的值
    # data_group['group'] = data_group['group'].replace(replacement_dict)


    if test_size_tmp != None:
        X_train, X_test, y_train, y_test = train_test_split(data_t, data_group, test_size=test_size_tmp, random_state=42)
    # %%

    # LGBM
    if model == "lgbm":
        lgbm = LGBMClassifier(objective='multiclass')
        lgbm.fit(X_train, y_train)
        y_pred = lgbm.predict(X_test)

        result1 = classification_report(y_test, y_pred,output_dict=True)
        print("LGBM Report:", )
        print(result1)
        result2 = accuracy_score(y_test, y_pred)
        print("Accuracy:", result2)
        
        y_pred_prob = lgbm.predict_proba(X_test)[:, 1]
        fpr, tpr, thresholds = roc_curve(y_test, y_pred_prob)
        roc_auc = auc(fpr, tpr)
        roc_data = pd.DataFrame({
            'FPR': fpr,
            'TPR': tpr,
            'Thresholds': thresholds
        })
        model_name = "LGBM"

    # randomForest
    elif model == "rf":
        rf = RandomForestClassifier(n_estimators=100)
        rf.fit(X_train, y_train)
        y_pred = rf.predict(X_test)

        result1 = classification_report(y_test, y_pred,output_dict=True)
        print("RandomForest Report:", )
        print(result1)
        result2 = accuracy_score(y_test, y_pred)
        print("Accuracy:", result2)
        
        y_pred_prob = rf.predict_proba(X_test)[:, 1]
        fpr, tpr, thresholds = roc_curve(y_test, y_pred_prob)
        roc_auc = auc(fpr, tpr)
        roc_data = pd.DataFrame({
            'FPR': fpr,
            'TPR': tpr,
            'Thresholds': thresholds
        })
        model_name = "RandomForest"
        
    # catBoost
    elif model == "catboost":
        cat_model = catboost.CatBoostClassifier(task_type="GPU",learning_rate=0.01)
        cat_model.fit(X_train, y_train)
        y_pred = cat_model.predict(X_test)

        result1 = classification_report(y_test, y_pred,output_dict=True)
        print("CatBoost Report:", )
        print(result1)
        result2 = accuracy_score(y_test, y_pred)
        print("Accuracy:", result2)

        y_pred_prob = cat_model.predict_proba(X_test)[:, 1]
        fpr, tpr, thresholds = roc_curve(y_test, y_pred_prob)
        roc_auc = auc(fpr, tpr)
        roc_data = pd.DataFrame({
            'FPR': fpr,
            'TPR': tpr,
            'Thresholds': thresholds
        })
        model_name = "CatBoost"
        
    # logistic
    elif model == "logistic":
        lr = LogisticRegression(solver="liblinear", max_iter=300)
        lr.fit(X_train, y_train)
        y_pred = lr.predict(X_test)

        result1 = classification_report(y_test, y_pred,output_dict=True)
        print("Logistic Report:", )
        print(result1)
        result2 = accuracy_score(y_test, y_pred)
        print("Accuracy:", result2)
        
        y_pred_prob = lr.predict_proba(X_test)[:, 1]
        fpr, tpr, thresholds = roc_curve(y_test, y_pred_prob)
        roc_auc = auc(fpr, tpr)
        roc_data = pd.DataFrame({
            'FPR': fpr,
            'TPR': tpr,
            'Thresholds': thresholds,
            'AUC value': roc_auc
        })
        model_name = "Logistic"

    # SVC
    elif model == "svc":
        svm = SVC(kernel="linear")
        svm.fit(X_train, y_train)
        y_pred = svm.predict(X_test)

        result1 = classification_report(y_test, y_pred,output_dict=True)
        print("SVC Report:", )
        print(result1)
        result2 = accuracy_score(y_test, y_pred)
        print("Accuracy:", result2)

        y_pred_prob = svm.predict_proba(X_test)[:, 1]
        fpr, tpr, thresholds = roc_curve(y_test, y_pred_prob)
        roc_auc = auc(fpr, tpr)
        roc_data = pd.DataFrame({
            'FPR': fpr,
            'TPR': tpr,
            'Thresholds': thresholds
        })
        model_name = "SVC"
        
    # DecisionTree
    elif model == "dt":
        dt = DecisionTreeClassifier()
        dt.fit(X_train, y_train)
        y_pred = dt.predict(X_test)

        result1 = classification_report(y_test, y_pred,output_dict=True)
        print("DecisionTree Report:", )
        print(result1)
        result2 = accuracy_score(y_test, y_pred)
        print("Accuracy:", result2)

        y_pred_prob = dt.predict_proba(X_test)[:, 1]
        fpr, tpr, thresholds = roc_curve(y_test, y_pred_prob)
        roc_auc = auc(fpr, tpr)
        roc_data = pd.DataFrame({
            'FPR': fpr,
            'TPR': tpr,
            'Thresholds': thresholds
        })
        model_name = "DecisionTree"
        
    # GaussianNB
    elif model == "gub":
        gnb = GaussianNB()
        gnb.fit(X_train, y_train)
        y_pred = gnb.predict(X_test)

        result1 = classification_report(y_test, y_pred,output_dict=True)
        print("GaussianNB Report:", )
        print(result1)
        result2 = accuracy_score(y_test, y_pred)
        print("Accuracy:", result2)

        y_pred_prob = gnb.predict_proba(X_test)[:, 1]
        fpr, tpr, thresholds = roc_curve(y_test, y_pred_prob)
        roc_auc = auc(fpr, tpr)
        roc_data = pd.DataFrame({
            'FPR': fpr,
            'TPR': tpr,
            'Thresholds': thresholds
        })
        model_name = "GaussianNB"
        
    # MultinomialNB
    elif model == "mnb":
        mnb = MultinomialNB()
        mnb.fit(X_train, y_train)
        y_pred = mnb.predict(X_test)

        result1 = classification_report(y_test, y_pred,output_dict=True)
        print("MultinomialNB Report:", )
        print(result1)
        result2 = accuracy_score(y_test, y_pred)
        print("Accuracy:", result2)

        y_pred_prob = mnb.predict_proba(X_test)[:, 1]
        fpr, tpr, thresholds = roc_curve(y_test, y_pred_prob)
        roc_auc = auc(fpr, tpr)
        roc_data = pd.DataFrame({
            'FPR': fpr,
            'TPR': tpr,
            'Thresholds': thresholds
        })
        model_name = "MultinomialNB"
        
    # BernoulliNB
    elif model == "bnb":
        bnb = BernoulliNB()
        bnb.fit(X_train, y_train)
        y_pred = bnb.predict(X_test)

        result1 = classification_report(y_test, y_pred,output_dict=True)
        print("BernoulliNB Report:", )
        print(result1)
        result2 = accuracy_score(y_test, y_pred)
        print("Accuracy:", result2)

        y_pred_prob = bnb.predict_proba(X_test)[:, 1]
        fpr, tpr, thresholds = roc_curve(y_test, y_pred_prob)
        roc_auc = auc(fpr, tpr)
        roc_data = pd.DataFrame({
            'FPR': fpr,
            'TPR': tpr,
            'Thresholds': thresholds
        })
        model_name = "BernoulliNB"
        
        
    report_df = pd.DataFrame(result1).transpose()
    accuracy_df = pd.DataFrame({'Accuracy': [result2]})

    return report_df, accuracy_df, roc_data, roc_auc, model_name
