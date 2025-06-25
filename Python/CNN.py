import tensorflow as tf

from keras.layers import Conv2D, MaxPooling2D, Flatten, Dense

from keras.models import Sequential

from keras.utils import to_categorical



# Carregar os dados

(x_train, y_train), (x_test, y_test) = tf.keras.datasets.mnist.load_data()



# Pré-processamento

x_train = x_train.reshape(-1, 28, 28, 1).astype("float32") / 255.0

x_test  = x_test.reshape(-1, 28, 28, 1).astype("float32") / 255.0

y_train = to_categorical(y_train, 10)

y_test  = to_categorical(y_test, 10)



# Definir a CNN

model = Sequential([

    Conv2D(32, kernel_size=(3, 3), activation='sigmoid', input_shape=(28, 28, 1)),

    MaxPooling2D(pool_size=(2, 2)),

    Conv2D(64, kernel_size=(3, 3), activation='sigmoid'),

    MaxPooling2D(pool_size=(2, 2)),

    Flatten(),

    Dense(128, use_bias=True, activation='relu'),

    Dense(10, use_bias=True, activation='softmax')

])



# Compilar o modelo

model.compile(optimizer='adam',

              loss='categorical_crossentropy',

              metrics=['accuracy'])



# Treinar o modelo

model.fit(x_train, y_train,

          epochs=5,

          batch_size=128,

          validation_split=0.2)



# Avaliar no conjunto de teste

test_loss, test_accuracy = model.evaluate(x_test, y_test)

print(f"Loss: {test_loss:.4f}, Accuracy: {test_accuracy:.4f}")

