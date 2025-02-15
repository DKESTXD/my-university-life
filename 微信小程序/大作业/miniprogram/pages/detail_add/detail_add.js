// pages/detail_add/detail_add.js
const app=getApp();
const date = new Date();
const years = []
const months = []
const days = []
const bigMonth = [1, 3, 5, 7, 8, 10, 12]
var getYear = date.getFullYear()
var getMonth = date.getMonth()
var getDate = date.getDate()
for (let i = getYear; i <= getYear + 40; i++) {years.push(i);}
for (let i = 1; i <= 12; i++) {months.push(i);}
for (let i = 1; i <= 31; i++) {days.push(i);}
Page({

  /**
   * 页面的初始数据
   */
  data: {
    titleValue:'',/*日程标题的值*/
    star1:true,star2:false,star3:false,star4:false,star5:false,
    urgent_level:1,/*紧急程度 默认1*/
    urgent_statue:'以后的事以后再说',
    inputLength:0,/*日程正文的字数*/
    textValue:'',/*日程正文*/
    reminder:false,/*是否提醒对应值 是为true*/
    settingValue:"仅仅本次",/*循环设置对应值*/
    showModalView:false,
    showModalViewCancel:false,
    hour: "00",/*截止时间的小时数*/
    minute:"00",/*截止时间的分钟数*/
    value: [0, getMonth, getDate - 1, "00", "00"],/*打开时间选择的时候显示的值*/
    hours: ["00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23"],
    minutes: ["00","01","02","03","04","05","06","07","08","09","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36","37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54","55","56","57","58","59"],
    timeInput: '',
    propDate: false,
    year: getYear,
    month: getMonth+1,
    day: getDate,
    years: years,
    months: months,
    days: days,
    setTime:""
  },
  /**
   * 生命周期函数--监听页面加载
   */
  onLoad(options) {
    var timelet=`${date.getFullYear()}-${(date.getMonth()+1)<10?`0${date.getMonth()+1}`:date.getMonth()+1}-${date.getDate()<10?`0${date.getDate()}`:date.getDate()} ${date.getHours()<10?`0${date.getHours()}`:date.getHours()}:${date.getMinutes()<10?`0${date.getMinutes()}`:date.getMinutes()}:00`;
    console.log(timelet);
    console.log(wx.getSetting());

  },
  refresh:function(){
    this.setData({
      showModalViewCancel: false,
      titleValue:'',
      star1:true,star2:false,star3:false,star4:false, star5:false,
    urgent_level:1,
    urgent_statue:'以后的事以后再说',
    inputLength:0,
    textValue:'',
    reminder:false,
    settingValue:"仅仅本次",
    hour: "00",
    minute:"00",
    value: [0, getMonth, getDate - 1, "00", "00"],
    timeInput: '',
    propDate: false,
    year: getYear,
    month: getMonth+1,
    day: getDate,
    });
  },
  testto:function(){
    var year=this.data.year;
    var month=this.data.month;
    var day=this.data.day;
    var hour=this.data.hour;
    var minute=this.data.minute;
    var timeset=`${year}-${month<10?`0${month}`:month}-${day<10?`0${day}`:day}T${hour}:${minute}:00`
    var timenow=`${date.getFullYear()}-${(date.getMonth()+1)<10?`0${date.getMonth()+1}`:date.getMonth()+1}-${date.getDate()<10?`0${date.getDate()}`:date.getDate()}T${date.getHours()<10?`0${date.getHours()}`:date.getHours()}:${date.getMinutes()<10?`0${date.getMinutes()}`:date.getMinutes()}:${date.getSeconds()<10?`0${date.getSeconds()}`:date.getSeconds()}`
    if(timeset<=timenow&&this.data.timeInput!=''){
      wx.showModal({
        title: '截止时间似乎已经过去了哦',
        icon:'warn',
        showCancel:false
      })
      console.log(timeset);
      console.log(timenow);
    }else if(this.data.reminder==true&&this.data.timeInput==''){
      wx.showModal({
        title: '您设置了提醒，请设置截止时间哦',
        icon:'warn',
        showCancel:false
      })
    }
    else if(this.data.titleValue==""){
      wx.showModal({
        title: '标题不能为空',
        icon:'warn',
        showCancel:false
      })
    }
    else{
      const db=wx.cloud.database({
        env:app.globalData.cloudID
      });
      if(this.data.reminder==true){
        db.collection("reminder").add({
          data:{
            titleValue:this.data.titleValue,
            textValue:this.data.textValue,
            urgent_level:this.data.urgent_level,
            settingValue:this.data.settingValue,
            timeInput:this.data.timeInput,
            done:false
          }
        })
      }
    wx.navigateTo({
      url: `../testpage/testpage`
    })
    db.collection("schedule").add({
      data:{
        titleValue:this.data.titleValue,
        textValue:this.data.textValue,
        urgent_level:this.data.urgent_level,
        reminder:this.data.reminder,
        settingValue:this.data.settingValue,
        timeInput:this.data.timeInput,
        setTime:this.data.setTime
      },
      success:res=>{
        wx.showToast({
          title: '新日程创建成功',
        })
        console.log("新增记录成功",res._id)
      },
      fail:err=>{
        wx.showToast({
          icon:noneParamsEaseFuncs,
          title: '新日程创建失败',
        })
        console.log("新增记录失败",err)
      }
    })
    this.refresh();
  }
  },
  titleInput:function(e){
    this.setData({
      titleValue: e.detail.value
    })
  },
  textInput: function(e) {
    this.setData({
      textValue:e.detail.value,
      inputLength: e.detail.value.length
    });
    if (e.detail.value.length > 300) {
      wx.showToast({
        title: '输入字数超出限制',
        icon: 'none',
        duration: 2000
      });
      this.setData({
        inputLength: 300
      });
    }
  },
  star1_tap:function(){
    this.setData({
      star1:true,star2:false,star3:false,star4:false,star5:false,urgent_statue:'以后的事以后再说',
      urgent_level:1
    });
  },
  star2_tap:function(){
    this.setData({
      star1:true,star2:true,star3:false,star4:false,star5:false,urgent_statue:'轻轻松 时间还很充足',
      urgent_level:2
    });
  },
  star3_tap:function(){
    this.setData({
      star1:true,star2:true,star3:true,star4:false,star5:false,urgent_statue:'你意识到你不能再拖了',
      urgent_level:3
    });
  },
  star4_tap:function(){
    this.setData({
      star1:true,star2:true,star3:true,star4:true,star5:false,urgent_statue:'非常重要 但是似乎会有比它更重要的事',urgent_level:4
    });
  },
  star5_tap:function(){
    this.setData({
      star1:true,star2:true,star3:true,star4:true,star5:true,urgent_statue:'🔥🔥🔥十万火急🔥🔥🔥',
      urgent_level:5
    });
  },
  reminder_change: function(e) {
    const value1 = e.detail.value; 
    const that = this; 
    if (value1 === true) {
      this.setData({
        reminder: true
    });
        wx.requestSubscribeMessage({
            tmplIds: ['v3J9x3kz6FYW3Xt4XwJ-d3qpbTU-J19XBh91HwOKqLA'],
            success(res) {
                if (res['v3J9x3kz6FYW3Xt4XwJ-d3qpbTU-J19XBh91HwOKqLA'] === 'accept') {
                    wx.showToast({
                        title: '成功订阅日程提醒',
                        icon: "none"
                    });
                } else if (res['v3J9x3kz6FYW3Xt4XwJ-d3qpbTU-J19XBh91HwOKqLA'] === 'reject') {
                    wx.showToast({
                        title: '您将无法使用日程提醒服务',
                        icon: "none"
                    });
                    that.setData({
                        reminder: false 
                    });
                }
            },
            fail(err) {
                console.error('订阅消息调用失败', err);
                wx.showToast({
                    title: '订阅失败，请检查设置',
                    icon: "none"
                });
                that.setData({
                    reminder: false 
                });
            }
        });
    } else {
        this.setData({
            reminder: false
        });
    }
},
  showModal: function() {
    this.setData({
      showModalView: true
    });
  },
  selectSetting_only: function(e) {
    this.setData({
      settingValue: "仅仅本次",
      showModalView:false
    });
  },
  selectSetting_weekly: function(e) {
    this.setData({
      settingValue: "每周一次",
      showModalView:false
    });
  },
  selectSetting_monthly: function(e) {
    this.setData({
      settingValue: "每月一次",
      showModalView:false
    });
  },
  cancel: function() {
    this.setData({
      showModalViewCancel: false
    });
  },
  cancel_tap: function() {
    this.setData({
      showModalViewCancel: true
    });
  },

  confirm_edit:function(){
    wx.navigateTo({
      url: `../testpage/testpage`
    });
    this.setData({
      showModalViewCancel: false,
    });
    this.refresh();
  },
  closePick() {
    this.setData({
      propDate: false
    })
  },
  openPick () {
    let valueTime = ''
    if(this.data.timeInput){
      const res = this.data.timeInput.split("-")
      const year = res[0].substring(2)
      const newRes = [ year-25, res[1]-1, res[2]-1,res[3],res[4]]
      valueTime = newRes
    }else{
      const year = this.data.year.toString().substring(2)
      valueTime = [year-25 , this.data.month-1 ,this.data.day-1,"00","00"]
    }
    this.setData({
      propDate: true,
      value:valueTime,
    })
    console.log(valueTime);
  },
  clearPick() {
    this.setData({
      timeInput:'',
      value:'',
    })
  },
  
  contains(arr, obj) {
    var i = arr.length;
    while (i--) {
      if (arr[i] === obj) {
        return true;
      }
    }
    return false;
  },
  setDays(day) {
    const temp = [];
    for (let i = 1; i <= day; i++) {
      temp.push(i)
    }
    this.setData({
      days: temp,
    })
  },
  //选择滚动器改变触发事件
  bindChange (e) {
    const val = e.detail.value;
    const setYear = this.data.years[val[0]];
    const setMonth = this.data.months[val[1]];
    const setDay = this.data.days[val[2]]
    const sethour = this.data.hours[val[3]];
    const setminute = this.data.minutes[val[4]];
    this.setData({
      year: setYear,
      month: setMonth,
      day: setDay,
      hour: sethour,
      minute: setminute,
    })
    if (setMonth === 2) {
      this.setDays((setYear % 4 === 0 && setYear % 100 !== 0) ? 29 : 28)
    } else {
      this.setDays(this.contains(bigMonth, setMonth)? 31 : 30)
    }
    this.setData({
      timeInput:setYear+ '-' +setMonth+ '-' +setDay+ '-' +sethour + '-' + setminute 
    })
  },
  /**
   * 生命周期函数--监听页面初次渲染完成
   */
  onReady() {

  },

  /**
   * 生命周期函数--监听页面显示
   */
  onShow() {

  },

  /**
   * 生命周期函数--监听页面隐藏
   */
  onHide() {

  },

  /**
   * 生命周期函数--监听页面卸载
   */
  onUnload() {

  },

  /**
   * 页面相关事件处理函数--监听用户下拉动作
   */
  onPullDownRefresh() {

  },

  /**
   * 页面上拉触底事件的处理函数
   */
  onReachBottom() {

  },

  /**
   * 用户点击右上角分享
   */
  onShareAppMessage() {

  }
})