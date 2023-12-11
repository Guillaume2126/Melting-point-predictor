FROM python:3.10.6-buster
COPY API/fast.py /API/fast.py
COPY Preprocessing/Create_dataframe.py /Preprocessing/Create_dataframe.py
RUN mkdir -p /Model/
COPY Model/best_full_data_model.joblib /Model/
COPY requirements.txt /requirements.txt
RUN pip install --upgrade pip
RUN pip install -r requirements.txt
CMD uvicorn API.fast:app --host 0.0.0.0 --port $PORT
