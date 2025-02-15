// app.js
App({
  onLaunch() {
    // 展示本地存储能力
    const logs = wx.getStorageSync('logs') || []
    logs.unshift(Date.now())
    wx.setStorageSync('logs', logs)
    if(!wx.cloud){
      console.error("初始化失败")
    }else{
      wx.cloud.init({
        env:"cloud1-3gpruo2f69271d50",
        traceUser:true,
      })
      wx.cloud.callFunction({
        name:"helloCloud",
        data: {
          message: 'helloCloud',
          }
      }).then(res=>{
        const user=res.result.event.userInfo;
        this.globalData.userInfo=user;
      })
    }
    // 登录
    wx.login({
      success: res => {
        // 发送 res.code 到后台换取 openId, sessionKey, unionId
      }
    })
  },
  globalData: {
    userInfo: null,
    cloudID:"cloud1-3gpruo2f69271d50"
  }
})
