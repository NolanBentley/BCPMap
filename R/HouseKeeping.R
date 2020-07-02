newObj<-function(env,baseObjectList){
    ls(name = env)[!ls(name = env)%in%baseObjectList]
}

