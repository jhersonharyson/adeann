import matplotlib.pyplot as plt
import pandas as pd
from keras.layers import Dense
from keras.models import Sequential
from sklearn.model_selection import train_test_split

limiar = 0.1

data = pd.read_csv('dados1.csv', encoding='utf-8', delimiter=';')
data_y = (data['Pesado'] + data['Medio'] + data['Leve'])/3
data_x = data.drop(["Pesado", "Medio", "Leve"], axis=1)
data_x = data.drop(["DiaInicio", "DiaFinal", "MesInicio", "AnoInicio", "MesFinal", "AnoFinal"], axis=1)

x_train, x_test, y_train, y_test = train_test_split(data_x, data_y, test_size=0.33,  random_state=42)

model = Sequential()
model.add(Dense(units=4, activation='elu', input_dim=8))
model.add(Dense(units=4, activation='elu'))
model.add(Dense(units=1, activation='linear'))

model.compile(loss='mean_absolute_error',
              optimizer='adam',
              metrics=['mean_absolute_error'])

model.fit(x_train, y_train, epochs=5000, batch_size=40)


predictions = model.predict(x_test)
total = len(x_test)
acertos = 0
rounded = [x for x in predictions] # [np.round(x) for x in predictions]
err = 0.0
for i, predic in enumerate(rounded):
    err = 0.5 * (y_test.iloc[i] - predic) ** 2
    print("\nPadrao>>", i)
    print("\ncalculado>>" + str(predic) + "   \tdesejado>>" + str(y_test.iloc[i]) + "  \tErro>>", err)

    if predic <= (y_test.iloc[i]  + y_test.iloc[i] *limiar) and predic >= (y_test.iloc[i]  - y_test.iloc[i] *limiar):
        acertos += 1


EMQ = err
y_wrapper = y_test.reset_index(drop=True)

FITNESS = (acertos / total)

print("\nemq>> ", EMQ)
print("\nfitness>> ", FITNESS)
print("\n\n<<Pesos Camada Entrada Oculta>>\n")

plt.plot(predictions, color="r")
plt.plot(y_wrapper, color="g")
plt.show()

# TODO: melhorar o resultado abaixo -------------------------------------------------------------------
#print(str(model.get_weights()))



