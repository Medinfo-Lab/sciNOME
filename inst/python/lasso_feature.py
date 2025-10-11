from sklearn.linear_model import Lasso
import pandas as pd
from sklearn.model_selection import train_test_split, GridSearchCV
import numpy as np




def lasso_feature(data_tmp,sample_data_tmp,test_size_tmp):
    data = data_tmp
    sample_data = sample_data_tmp

    #%%
    data_t = data.T
    # DEG_data_t_tmp = pd.concat([DEG_data_t,DEG_sample_data],axis=1)
    # DEG_sample_data_encoded = pd.get_dummies(DEG_sample_data, columns=['group'])
    # df_encoded = pd.get_dummies(DEG_data_t_tmp, columns=['group'])  # 对'Tumor'列进行独热编码

    # 创建一个映射字典，将'Normal'映射为0，'Tumor'映射为1
    # replacement_dict = {'normal': 0, 'tumor': 1}

    # 使用.replace()方法替换'Sample'列中的值
    # sample_data['group_id'] = sample_data['group'].replace(replacement_dict)
    # sample_data['group'] = sample_data['group'].replace(group_mapping)

    X_train, X_test, y_train, y_test = train_test_split(data_t, sample_data, test_size = test_size_tmp, random_state=42)

  
    #%%
    # 设置LASSO回归模型
    lasso = Lasso()

    # 使用网格搜索找到最佳的正则化参数alpha
    param_grid = {'alpha': np.logspace(-4, 4, 50)}
    lasso_cv = GridSearchCV(lasso, param_grid, cv=5)
    lasso_cv.fit(X_train, y_train)

    print("Best alpha:", lasso_cv.best_params_['alpha'])

    # 使用最佳的正则化参数重新训练LASSO模型
    lasso_best = Lasso(alpha=lasso_cv.best_params_['alpha'])
    lasso_best.fit(X_train, y_train)

    # 获取特征的系数
    coef = lasso_best.coef_

    # 选择系数非零的特征
    selected_features = np.where(coef != 0)[0]
    selected_genes = data_t.columns[selected_features]

    # 输出选定的特征基因
    # print("Selected genes:", selected_genes.tolist())


    #%%
    selected_genes_data = pd.DataFrame(selected_genes)
    data_t_model = data_t[selected_genes_data.iloc[:,0]]
    return data_t_model
